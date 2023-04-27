class Error:
    error_code = 0
    error_text = ""
    error_name = "Error"
    
    def __init__(self, error_code: int, error_text: str, error_name: str = ""):
        self.error_code = error_code
        self.error_text = error_text
    
    def get_error_code(self):
        return self.error_code
    
    def get_error_text(self):
        return self.error_text
    
    def get_error_name(self):
       return self.error_name 
    
class StatisticsError(Error):
    base_error_code = 3000
    error_code = 3000
    error_name = "StatisticsError"
    error_text = ""
    
    """
    A StatisticsError, used with fischer and other statistical computations during the scoring of products.
    Base error code for StatisticsError is 3000, you can add 'error_code_specifiers' (eg. 1, 2, 3) to create more specific error codes (self.error_code = self.base_error_code + error_code_specific)
    """
    def __init__(self, error_text: str, error_code_specific: int = 0, error_name: str = ""):
        self.error_text = error_text
        self.error_code = self.base_error_code + error_code_specific
        self.error_name = error_name if error_name != "" else self.error_name
    