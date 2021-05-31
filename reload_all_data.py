from bookkeeping import bookkeeping


def reload_all_data():

    reloader = bookkeeping()
    reloader.update_from_scratch()
    b = bookkeeping()
    b.make_cls_crossed_data()
# reload_all_data()
