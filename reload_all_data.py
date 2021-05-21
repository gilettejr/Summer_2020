from bookkeeping import bookkeeping
def reload_all_data():
    
    reloader=bookkeeping()
    reloader.update_from_scratch()
#reload_all_data()