import os

def simplify_path(path,remove_gz=True):
    """Removes dir and extension from a filepath.
        checks if file has an e
    """
    name,ext= os.path.splitext(os.path.basename(path))

    if remove_gz & (ext=='.gz'):
        name=  os.path.splitext(name)[0]

    return name
