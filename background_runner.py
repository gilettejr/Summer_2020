from background_constructor import background_constructor
from data_processor import data_processor
from close_data_processor import close_data_processor

class background_runner:
    
    def __init__(self):

        def make_close_backgrounds(galaxy,stars):
            
            galaxy_object=background_constructor(galaxy=galaxy,stage='cm')
            galaxy_object.find_background_grad(stars=stars)
            galaxy_object.fit_close_background()
            galaxy_object.construct_slices(stars=stars)
            galaxy_object.find_close_slice_profile()
            
        self.make_close_backgrounds=make_close_backgrounds
    
    def make_all_close_backgrounds(self):
        
        make_close_backgrounds=self.make_close_backgrounds
        
        star_list=['agb','m','c']
        
        galaxy_list=['ngc205','m32']
        
        for i in galaxy_list:
            
            for j in star_list:
                
                make_close_backgrounds(galaxy=i,stars=j)
    
