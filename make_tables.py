from data_readall import data_readall
from data_load import data_loader
import pandas as pd
import matplotlib.pyplot as plt


def make_table():
    processed_stages = ['cls_cut', 'fore_cut', 'cm']
    # galaxies=['ngc147','ngc185','ngc205','m32']

    dr = data_readall(stage='cls_cut')
    n147_cls = dr.n147.data
    n185_cls = dr.n185.data
    n205_cls = dr.n205.data
    m32_cls = dr.m32.data

    dr = data_readall(stage='cls_crossed')
    n147_cls_crossed = dr.n147.data
    n185_cls_crossed = dr.n185.data
    n205_cls_crossed = dr.n205.data
    m32_cls_crossed = dr.m32.data

    dr = data_readall(stage='fore_cut')
    n147_fore = dr.n147.data
    n185_fore = dr.n185.data
    n205_fore = dr.n205.data
    m32_fore = dr.m32.data

    dr = data_readall(stage='cm')
    n147_c = dr.n147.cdata
    n185_c = dr.n185.cdata
    n205_c = dr.n205.cdata
    m32_c = dr.m32.data

    n147_m = dr.n147.mdata
    n185_m = dr.n185.mdata
    n205_m = dr.n205.mdata
    m32_m = dr.m32.mdata

    n147_agb = dr.n147.data
    n185_agb = dr.n185.data
    n205_agb = dr.n205.data
    m32_agb = dr.m32.data

    n147_fore = pd.concat([n147_fore, n147_cls]).drop_duplicates(keep=False)
    n185_fore = pd.concat([n185_fore, n185_cls]).drop_duplicates(keep=False)
    n205_fore = pd.concat([n205_fore, n205_cls]).drop_duplicates(keep=False)
    m32_fore = pd.concat([m32_fore, m32_cls]).drop_duplicates(keep=False)

    n147_rgb = pd.concat([n147_fore, n147_agb]).drop_duplicates(keep=False)
    n185_rgb = pd.concat([n185_fore, n185_agb]).drop_duplicates(keep=False)
    n205_rgb = pd.concat([n205_fore, n205_agb]).drop_duplicates(keep=False)
    m32_rgb = pd.concat([m32_fore, m32_agb]).drop_duplicates(keep=False)

    n147_gai = pd.concat([n147_cls, n147_cls_crossed]
                         ).drop_duplicates(keep=False)
    n185_gai = pd.concat([n185_cls, n185_cls_crossed]
                         ).drop_duplicates(keep=False)
    n205_gai = pd.concat([n205_cls, n205_cls_crossed]
                         ).drop_duplicates(keep=False)
    m32_gai = pd.concat([m32_cls, m32_cls_crossed]
                        ).drop_duplicates(keep=False)

    n147_gai = pd.concat([n147_gai, n147_fore]
                         ).drop_duplicates(keep=False)
    n185_gai = pd.concat([n185_gai, n185_fore]
                         ).drop_duplicates(keep=False)
    n205_gai = pd.concat([n205_gai, n205_fore]
                         ).drop_duplicates(keep=False)
    m32_gai = pd.concat([m32_gai, m32_fore]
                        ).drop_duplicates(keep=False)

    n147_gai = pd.concat([n147_gai, n147_rgb]
                         ).drop_duplicates(keep=False)
    n185_gai = pd.concat([n185_gai, n185_rgb]
                         ).drop_duplicates(keep=False)
    n205_gai = pd.concat([n205_gai, n205_rgb]
                         ).drop_duplicates(keep=False)
    m32_gai = pd.concat([m32_gai, m32_rgb]
                        ).drop_duplicates(keep=False)

    dl = data_loader(galaxy='ngc147', CLS=False, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    n147_raw = dl.data
    dl = data_loader(galaxy='ngc185', CLS=False, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    n185_raw = dl.data
    dl = data_loader(galaxy='ngc205', CLS=False, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    n205_raw = dl.data
    dl = data_loader(galaxy='m32', CLS=False, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    m32_raw = dl.data
    dl = data_loader(galaxy='ngc147', CLS=True, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    n147_mag = dl.data
    dl = data_loader(galaxy='ngc185', CLS=True, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    n185_mag = dl.data
    dl = data_loader(galaxy='ngc205', CLS=True, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    n205_mag = dl.data
    dl = data_loader(galaxy='m32', CLS=True, cls_bands='norm',
                     mag=False, ext=True, path_to_file='initial_data/')
    m32_mag = dl.data

    n147_mag_fail = pd.concat([n147_cls, n147_mag]
                              ).drop_duplicates(keep=False)
    n185_mag_fail = pd.concat([n185_cls, n185_mag]
                              ).drop_duplicates(keep=False)
    n205_mag_fail = pd.concat([n205_cls, n205_mag]
                              ).drop_duplicates(keep=False)
    m32_mag_fail = pd.concat([m32_cls, m32_mag]
                             ).drop_duplicates(keep=False)

    n147_cls_fail = pd.concat([n147_mag, n147_raw]
                              ).drop_duplicates(keep=False)
    n185_cls_fail = pd.concat([n185_mag, n185_raw]
                              ).drop_duplicates(keep=False)
    n205_cls_fail = pd.concat([n205_mag, n205_raw]
                              ).drop_duplicates(keep=False)
    m32_cls_fail = pd.concat([m32_mag, m32_raw]
                             ).drop_duplicates(keep=False)

    def add_flag(frame, flag, galaxy):
        frame['class'] = flag
        frame['galaxy'] = galaxy
        frame = frame.dropna()

    add_flag(n147_c, 'C-AGB', 'NGC147')
    add_flag(n185_c, 'C-AGB', 'NGC185')
    add_flag(n205_c, 'C-AGB', 'NGC205')
    add_flag(m32_c, 'C-AGB', 'M32')

    add_flag(n147_m, 'M-AGB', 'NGC147')
    add_flag(n185_m, 'M-AGB', 'NGC185')
    add_flag(n205_m, 'M-AGB', 'NGC205')
    add_flag(m32_m, 'M-AGB', 'M32')

    add_flag(n147_rgb, 'RGB', 'NGC147')
    add_flag(n185_rgb, 'RGB', 'NGC185')
    add_flag(n205_rgb, 'RGB', 'NGC205')
    add_flag(m32_rgb, 'RGB', 'M32')

    add_flag(n147_rgb, 'RGB', 'NGC147')
    add_flag(n185_rgb, 'RGB', 'NGC185')
    add_flag(n205_rgb, 'RGB', 'NGC205')
    add_flag(m32_rgb, 'RGB', 'M32')

    add_flag(n147_gai, 'FORE_GAIA', 'NGC147')
    add_flag(n185_gai, 'FORE_GAIA', 'NGC185')
    add_flag(n205_gai, 'FORE_GAIA', 'NGC205')
    add_flag(m32_gai, 'FORE_GAIA', 'M32')

    add_flag(n147_fore, 'FORE_SEQ', 'NGC147')
    add_flag(n185_fore, 'FORE_SEQ', 'NGC185')
    add_flag(n205_fore, 'FORE_SEQ', 'NGC205')
    add_flag(m32_fore, 'FORE_SEQ', 'M32')

    add_flag(n147_cls_fail, 'NOISE-CLS', 'NGC147')
    add_flag(n185_cls_fail, 'NOISE-CLS', 'NGC185')
    add_flag(n205_cls_fail, 'NOISE-CLS', 'NGC205')
    add_flag(m32_cls_fail, 'NOISE-CLS', 'M32')

    add_flag(n147_mag_fail, 'NOISE-MAG', 'NGC147')
    add_flag(n185_mag_fail, 'NOISE-MAG', 'NGC185')
    add_flag(n205_mag_fail, 'NOISE-MAG', 'NGC205')
    add_flag(m32_mag_fail, 'NOISE-MAG', 'M32')

    table = pd.concat([n147_cls_fail, n185_cls_fail,
                       n205_cls_fail, m32_cls_fail, n147_mag_fail, n185_mag_fail,
                       n205_mag_fail, m32_mag_fail, n147_fore, n185_fore,
                       n205_fore, m32_fore, n147_rgb, n185_rgb,
                       n205_rgb, m32_rgb, n147_gai, n185_gai,
                       n205_gai, m32_gai, n147_m, n185_m,
                       n205_m, m32_m, n147_c, n185_c,
                       n205_c, m32_c])

    table.to_parquet('master_table')


make_table()
