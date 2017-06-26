from os.path import expanduser
home = expanduser("~")
import logging
from lockfile import LockFile
from lockfile import LockTimeout
import numpy as np
import time
import unittest
import p_mymodule as my
logging.getLogger().setLevel(logging.DEBUG) #change the level of the display.

class locker_queue:
    locker_list=None
    cur_time=time.strftime("%Y%m%d_%H%M%S")

    def __init__(self):
        self.locker_list=[]

    def get_new_locker(self, locker_name):
        new_locker_file=self.cur_time +"." + locker_name + "." + my.f_id_generator(5)
        self.locker_list.append(new_locker_file)
        return new_locker_file

    def check_lokcer(self, pattern, timeout=10000000):
    
        import os.path
        import time
        
        #import ipdb; ipdb.set_trace()
        enquiry_lockers = my.grep_list(pattern, self.locker_list)
        return_status = True
        
        if len(enquiry_lockers) == 0:
            logging.warning("The equiery lokcer doesn't exist :" + pattern)
            return_status = False
        else:
            for locker_file in enquiry_lockers:
                print "Check the locker:" + locker_file
                total_time = 0
                while True:
                    if total_time > timeout:
                        return_status = False
                        logging.warning(locker_file + " reached timeout")
                        break
                    
                    if os.path.isfile("/homed/home/shi/locker_dir/%s"%(locker_file)):
                        print "Get the locker"
                        self.locker_list.remove(locker_file)
                        break
    
                    time.sleep(10)
                    total_time =  total_time + 10
            

        return return_status

class TestLockerFunctions(unittest.TestCase):
    def setUp(self):
        self.locker_list =  locker_queue()
        self.server_name = my.f_get_server_name()
        
    def test_add_locker(self):
        new_locker_file=self.locker_list.get_new_locker("abc")
        self.assertEqual([new_locker_file], self.locker_list.locker_list)

    def test_check_locker(self):
        locker1=self.locker_list.get_new_locker("aaa")
        locker2=self.locker_list.get_new_locker("bbb")
        locker3=self.locker_list.get_new_locker("ccc")
        my.f_write_locer_file(locker1,self.server_name)
        my.f_write_locer_file(locker2,self.server_name)
        logging.debug(self.locker_list.locker_list)
        self.assertTrue( self.locker_list.check_lokcer(".*(aaa|bbb)") )
        self.assertTrue( self.locker_list.check_lokcer(".*ddd") == False )
        self.assertTrue( self.locker_list.check_lokcer(".*ccc", timeout=3) == False)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase( TestLockerFunctions )
    unittest.TextTestRunner(verbosity=1,stream=sys.stderr).run( suite )

