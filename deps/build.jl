using PyCall
using Conda

try
    py"""
    from pip import main

    main(['install', '-y', 'matplotlib'])
    """
catch
    Conda.add("matplotlib")
end
