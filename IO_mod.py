# -*- coding: utf-8 -*-



# -----------------------------------------------------------------------------
def error(message):
    full_message = 'ERROR: ' + message
#    sys.exit(full_message)
    raise ValueError(full_message)
