from close_data_processor import close_data_processor
from background_constructor import background_constructor
from data_processor import data_processor


class FEH_finders:
    @staticmethod
    def find_far_FEH(galaxy):

        galaxy_object = data_processor(galaxy=galaxy)
        galaxy_object.construct_slices()
        background_object = background_constructor(galaxy=galaxy)
        m_background = background_object.find_background_density_border(
            stars='m')
        c_background = background_object.find_background_density_border(
            stars='c')
        galaxy_object.find_FEH_slices(m_background, c_background)
        galaxy_object.plot_ellipses(ellipticity=galaxy_object.ellipticity, a_inner=galaxy_object.a_width /
                                    2, a_outer=galaxy_object.outer_rad, PA=galaxy_object.rotation)

    @staticmethod
    def find_close_FEH(galaxy):

        print('You need to have run make_all_close_backgrounds in the background_runners class for this to work!!!')

        galaxy_object = close_data_processor(galaxy=galaxy)
        galaxy_object.construct_slices(galaxy=galaxy)
        galaxy_object.get_close_FEH_slices()
        galaxy_object.plot_ellipses(ellipticity=galaxy_object.ellipticity, a_inner=galaxy_object.a_width /
                                    2, a_outer=galaxy_object.outer_rad, PA=galaxy_object.rotation)
