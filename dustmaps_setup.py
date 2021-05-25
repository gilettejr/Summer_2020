import dustmaps.sfd
from dustmaps.config import config
# insert your preferred file path for map storage here
config['data_dir'] = '/home/robertsoncl/ifa/'

dustmaps.sfd.fetch()
