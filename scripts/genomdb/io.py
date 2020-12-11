import os
import gzip as gz


def simplify_path(path, remove_gz=True):
    """Removes dir and extension from a filepath.
        checks if file has an e
    """
    name, ext = os.path.splitext(os.path.basename(path))

    if remove_gz & (ext == ".gz"):
        name = os.path.splitext(name)[0]

    return name


def simply_open(filename, mode="r", *args, **kwargs):
    """open file irrespective if gz compressed or not"""

    if filename.endswith(".gz"):

        # To read file in textmode
        if mode in ["r", "a", "w", "x"]:
            mode += "t"

        return gz.open(filename, mode, *args, **kwargs)
    else:
        return open(filename, mode, *args, **kwargs)


def cat_files(files, outfilename, gzip=False):
    """ cat files in python
    """
    import shutil

    if gzip:
        outhandle = gz.open
    else:
        outhandle = open

    with outhandle(outfilename, "wb") as f_out:
        for f in files:
            with open(f, "rb") as f_in:
                shutil.copyfileobj(f_in, f_out)



def symlink_relative(files, input_dir, output_dir):
    """create symlink with and adjust for relative path"""

    input_dir_rel = os.path.relpath(input_dir, output_dir)

    for f in files:
        os.symlink(os.path.join(input_dir_rel, f), os.path.join(output_dir, f))
