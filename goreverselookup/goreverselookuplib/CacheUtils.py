from .Metrics import Metrics,basic_mirna_score
import json
import os
import logging
import time
from .JsonUtil import JsonUtil
from .Timer import Timer
import asyncio
import aiohttp

logger = logging.getLogger(__name__)

class Cacher():
    CACHE_FILEPATH_URLS = "" # filepath to the file containing online url queries and the URL RESPONSES
    CACHE_FILEPATH_UNIPROT = "" # filepath to the file containing uniprot api queries and their final results (after processing of the url responses)
    CACHE_FILEPATH_ENSEMBL = "" # filepath to the file containing ensembl api queries and their final results (after processing of the url responses)
    cached_urls = {}
    cached_uniprot = {}
    cached_ensembl = {}
    
    @classmethod
    def init(cls):
        """
        Initialises ConnectionCacher. This function must be called at the program startup in order to read
        old urls into the cls.cached_urls dictionary.

        Usage:
            model = ReverseLookup.load_model("diabetes_angio_4/model_async_test.json") # make sure that model products are already computed
            Cacher.init()
            fetch_ortholog_products(refetch=True, run_async=False)
        """
        cls.CACHE_FILEPATH_URLS = "cache/connection_cache.json"
        cls.CACHE_FILEPATH_UNIPROT = "cache/uniprot_cache.json"
        cls.CACHE_FILEPATH_ENSEMBL = "cache/ensembl_cache.json"
        cls.cached_urls = JsonUtil.load_json(cls.CACHE_FILEPATH_URLS)
        cls.cached_uniprot = JsonUtil.load_json(cls.CACHE_FILEPATH_UNIPROT)
        cls.cached_ensembl = JsonUtil.load_json(cls.CACHE_FILEPATH_ENSEMBL)

    @classmethod
    def store_data(cls, data_location:str, data_key:str, data_value, timestamp:str=""):
        """
        Stores 
            {
            "data_key": 
                "data_value": data_value,
                "timestamp": timestamp
            } 
        
        inside a particular json cache file, based on 'data_location'.
        The data_location options are:
          - "url": -> filepath = cache/connection_cache.json
          - "uniprot" -> filepath = cache/uniprot_cache.json
          - "ensembl" -> filepath = cache/ensembl_cache.json
        
        Url storage is intended for intermediate url query responses. Consider the following url query: f"https://rest.ensembl.org/homology/symbol/{species}/{id_url}?target_species=human;type=orthologues;sequence=none":
        Without request caching, this is the code algorithm:

            url = ...
            response = (Session).get(url)
            response_json = response.json()
        
        With request caching, the code algorithm is slightly modified:

            url = ...
            previous_response = Cacher.get_data(data_location="urls", data_key=url)
            if previous_response != None:
                response_json = previous_response
            else:
                response = (Session).get(url)
                response_json = response.json()
                Cacher.store_data(data_location="urls", data_key=url, data_value=response_json)
            # process response_json
        
        Alternatively, in the case of storage the response jsons of the queried urls, you can use ConnectionCacher:

            url = ...
            previous_response = ConnectionCacher.get_url_response(url)
            if previous_response != None:
                response_json = previous_response
            else:
                response = (Session).get(url)
                response_json = response.json()
                ConnectionCacher.store_url(url, response=response_json)
            # process response_json
        
        With this code, if the algorithm encounters and already queried url, it will pull its old response,
        rather than query a new one.
        """
        
        cached_data = {}
        # determine type of cached data
        match data_location:
            case "url":
                cached_data = cls.cached_urls
            case "uniprot":
                cached_data = cls.cached_uniprot
            case "ensembl":
                cached_data = cls.cached_ensembl

        # calculate current time
        if timestamp == "":
            timestamp = Timer.get_current_time()
        
        # update cached_data
        if data_key not in cached_data:
            cached_data[data_key] = {"data_value": data_value, "timestamp": timestamp}
        else: # this data_key already exists in previous data
            previous_data_timestamp = cached_data[data_key]["timestamp"]
            if Timer.compare_time(previous_data_timestamp, timestamp) == True: # will return true, if timestamp > previous_url_timestamp (timestamp is logged later in time than previous_url_timestamp)
                if data_value != None:
                    cached_data[data_key] = {"data_value": data_value, "timestamp": timestamp}
        
        # update class values with new cached data and save
        match data_location:
            case "url":
                cls.cached_urls = cached_data
                JsonUtil.save_json(cls.cached_urls, cls.CACHE_FILEPATH_URLS)
            case "uniprot":
                cls.cached_uniprot = cached_data
                JsonUtil.save_json(cls.cached_uniprot, cls.CACHE_FILEPATH_UNIPROT)
            case "ensembl":
                cls.cached_ensembl = cached_data
                JsonUtil.save_json(cls.cached_ensembl, cls.CACHE_FILEPATH_ENSEMBL)
    
    @classmethod
    def get_data(cls, data_location:str, data_key:str):
        cached_data = {}
        # determine type of cached data
        match data_location:
            case "url":
                cached_data = cls.cached_urls
            case "uniprot":
                cached_data = cls.cached_uniprot
            case "ensembl":
                cached_data = cls.cached_ensembl
        
        if cached_data != {}:
            if data_key in cached_data:
                logger.info(f"Successfully cached old data for {data_key}.")
                return cached_data[data_key]["data_value"]
            else:
                return None
        else:
            return None


class ConnectionCacher(Cacher):
    """
    ConnectionCacher accesses root/cache/connection_cache.json in order to store and retrieve old
    url connections and their belonging data (url response, time of request)

    A newer implementation, which combines connection caching, as well as caching of processed uniprot or
    ensembl function results, is the Cacher class. It is advisable to use the Cacher class with the
    parameter "url" as data_location in place of ConnectionCacher.

    NOTE: We could make three implementations of Cacher -> ConnectionCacher, UniprotCacher, EnsemblCacher,
    but that would add too complex functionality, which can be reasonably implemented in a single class.
    """
    CACHE_FILEPATH = "cache/connection_cache.json"
    cached_urls = {}

    @classmethod
    def init(cls):
        """
        Initialises ConnectionCacher. This function must be called at the program startup in order to read
        old urls into the cls.cached_urls dictionary.
        """
        cls.CACHE_FILEPATH = "cache/connection_cache.json"
        cls.cached_urls = JsonUtil.load_json(cls.CACHE_FILEPATH)

    @classmethod
    def store_url(cls, url:str, response, timestamp:str=""):
        """
        Stores the 'url' as the key, it's value is a dictionary comprised of 'response' and 'timestamp'.
        The key-value pair is stored in root/cache/connection_cache.json. If timestamp is not provided, then
        a timestamp will be calculated inside this function.

        Json outline:
        {url1 -> {"response": response1, "timestamp": timestamp1}},
        {url2 -> {"response": response2, "timestamp": timestamp2}},
        ...
        """
        # data = JsonUtil.load_json(cls.CACHE_FILEPATH)
        
        data = cls.cached_urls
        if timestamp == "":
            timestamp = Timer.get_current_time()

        if "url" not in data: # url doesn't exist in previous data -> add it
            data[url] = {"response": response, "timestamp": timestamp} # add new element
            cls.cached_urls = data # update cached urls
            JsonUtil.save_json(cls.cached_urls, cls.CACHE_FILEPATH) # save cached urls
        else: # this url already exists in previous data
            # previous_url_response = data[url]["response"]
            previous_url_timestamp = data[url]["timestamp"]
            if Timer.compare_time(previous_url_timestamp, timestamp) == True: # will return true, if timestamp > previous_url_timestamp (timestamp is logged later in time than previous_url_timestamp)
                if response != None:
                    data[url] = {"response": response, "timestamp": timestamp} # add new element
                    cls.cached_urls = data # update cached urls
                    JsonUtil.save_json(cls.cached_urls, cls.CACHE_FILEPATH) # save cached urls

    @classmethod
    def get_url_response(cls, url:str):
        """
        Obtains the response of the 'url' from previously cached urls, if the same url already exists.
        Previously cached urls and their responses are stored in root/cache/connection_cache.json.

        Returns None either if url doesn't exist or if the response of the url is stored as None.
        """
        if cls.cached_urls == {}:
            logger.warning(f"Cached_urls variable is empty! Did you forget to call ConnectionCacher.init()?")
            cls.init()
        
        if cls.cached_urls != {}:
            if url in cls.cached_urls:
                # TODO: implement url age option, eg. if the user selects "previous month", if the url is older than that, return None
                logger.info(f"Cached response for {url}: {cls.cached_urls[url]['response']}")
                return cls.cached_urls[url]["response"]
            else: # url wasn't found
                return None
        else:
            return None
            



    