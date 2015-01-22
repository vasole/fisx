import os
# this will be filled by the setup
FISX_DATA_DIR = 'DATA_DIR_FROM_SETUP'
# This is to be filled by the setup
FISX_DOC_DIR = 'DOC_DIR_FROM_SETUP'
# what follows is only used in frozen versions
if not os.path.exists(FISX_DATA_DIR):
    tmp_dir = os.path.dirname(__file__)
    old_tmp_dir = tmp_dir + "dummy"
    basename = os.path.basename(FISX_DATA_DIR)
    FISX_DATA_DIR = os.path.join(tmp_dir, "fisx", basename)
    while (len(FISX_DATA_DIR) > 14) and (tmp_dir != old_tmp_dir):
        if os.path.exists(FISX_DATA_DIR):
            break
        old_tmp_dir = tmp_dir
        tmp_dir = os.path.dirname(tmp_dir)
        FISX_DATA_DIR = os.path.join(tmp_dir, basename)

# this is used in build directory
if not os.path.exists(FISX_DATA_DIR):
    tmp_dir = os.path.dirname(__file__)
    old_tmp_dir = tmp_dir + "dummy"
    basename = os.path.basename(FISX_DATA_DIR)
    FISX_DATA_DIR = os.path.join(tmp_dir, "fisx", basename)
    while (len(FISX_DATA_DIR) > 14) and (tmp_dir != old_tmp_dir):
        if os.path.exists(FISX_DATA_DIR):
            break
        old_tmp_dir = tmp_dir
        tmp_dir = os.path.dirname(tmp_dir)
        FISX_DATA_DIR = os.path.join(tmp_dir, "fisx", basename)

if not os.path.exists(FISX_DATA_DIR):
    raise IOError('%s directory not found' % basename)

# do the same for the directory containing HTML files
if not os.path.exists(FISX_DOC_DIR):
    tmp_dir = os.path.dirname(__file__)
    old_tmp_dir = tmp_dir + "dummy"
    basename = os.path.basename(FISX_DOC_DIR)
    FISX_DOC_DIR = os.path.join(tmp_dir,basename)
    while (len(FISX_DOC_DIR) > 14) and (tmp_dir != old_tmp_dir):
        if os.path.exists(FISX_DOC_DIR):
            break
        old_tmp_dir = tmp_dir
        tmp_dir = os.path.dirname(tmp_dir)
        FISX_DOC_DIR = os.path.join(tmp_dir, basename)
if not os.path.exists(FISX_DOC_DIR):
    raise IOError('%s directory not found' % basename)

if not os.path.exists(FISX_DOC_DIR):
    tmp_dir = os.path.dirname(__file__)
    old_tmp_dir = tmp_dir + "dummy"
    basename = os.path.basename(FISX_DOC_DIR)
    FISX_DOC_DIR = os.path.join(tmp_dir,basename)
    while (len(FISX_DOC_DIR) > 14) and (tmp_dir != old_tmp_dir):
        if os.path.exists(FISX_DOC_DIR):
            break
        old_tmp_dir = tmp_dir
        tmp_dir = os.path.dirname(tmp_dir)
        FISX_DOC_DIR = os.path.join(tmp_dir, "fisx", basename)
if not os.path.exists(FISX_DOC_DIR):
    raise IOError('%s directory not found' % basename)
