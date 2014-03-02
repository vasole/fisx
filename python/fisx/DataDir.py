import os
# this will be filled by the setup
DATA_DIR = 'DATA_DIR_FROM_SETUP'
# This is to be filled by the setup
DOC_DIR = 'DOC_DIR_FROM_SETUP'
# what follows is only used in frozen versions
if not os.path.exists(DATA_DIR):
    tmp_dir = os.path.dirname(__file__)
    basename = os.path.basename(DATA_DIR)
    DATA_DIR = os.path.join(tmp_dir,basename)
    while len(DATA_DIR) > 14:
        if os.path.exists(DATA_DIR):
            break
        tmp_dir = os.path.dirname(tmp_dir)
        DATA_DIR = os.path.join(tmp_dir, basename)

# this is used in build directory
if not os.path.exists(DATA_DIR):
    tmp_dir = os.path.dirname(__file__)
    basename = os.path.basename(DATA_DIR)
    DATA_DIR = os.path.join(tmp_dir, "fisx", basename)
    while len(DATA_DIR) > 14:
        if os.path.exists(DATA_DIR):
            break
        tmp_dir = os.path.dirname(tmp_dir)
        DATA_DIR = os.path.join(tmp_dir, "fisx", basename)

if not os.path.exists(DATA_DIR):
    raise IOError('%s directory not found' % basename)
# do the same for the directory containing HTML files
if not os.path.exists(DOC_DIR):
    tmp_dir = os.path.dirname(__file__)
    basename = os.path.basename(DOC_DIR)
    DOC_DIR = os.path.join(tmp_dir,basename)
    while len(DOC_DIR) > 14:
        if os.path.exists(DOC_DIR):
            break
        tmp_dir = os.path.dirname(tmp_dir)
        DOC_DIR = os.path.join(tmp_dir, basename)
if not os.path.exists(DOC_DIR):
    raise IOError('%s directory not found' % basename)

if not os.path.exists(DOC_DIR):
    tmp_dir = os.path.dirname(__file__)
    basename = os.path.basename(DOC_DIR)
    DOC_DIR = os.path.join(tmp_dir, "fisx", basename)
    while len(DOC_DIR) > 14:
        if os.path.exists(DOC_DIR):
            break
        tmp_dir = os.path.dirname(tmp_dir)
        DOC_DIR = os.path.join(tmp_dir, "fisx", basename)
if not os.path.exists(DOC_DIR):
    raise IOError('%s directory not found' % basename)
