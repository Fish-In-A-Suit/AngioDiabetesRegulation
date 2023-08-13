# This test file demonstrates the use of ConnectionCacher.

from goreverselookuplib.CacheUtils import ConnectionCacher

test_url = "test_url.com"
test_url_response = "test response"
ConnectionCacher.init()
ConnectionCacher.store_url(test_url, test_url_response)
ConnectionCacher.get_url_response(test_url)