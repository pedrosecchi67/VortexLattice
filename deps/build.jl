using PyCall

py"""
try:
    from pip import main

    main(['install', '-y', 'matplotlib'])
except:
    print("Unable to install matplotlib. Plotting functions may be unavailable for execution, as will be warned during import")
"""
