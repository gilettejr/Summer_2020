from dustmaps.config import config
#insert your preferred file path for map storage here
config['data_dir'] = 'file/path/to/put/maps/in'

import dustmaps.sfd
dustmaps.sfd.fetch()