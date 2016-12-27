import os
import re


def ensure_requisite_folders(path):
    folder = os.path.split(path)[0]
    if len(folder) and not os.path.exists(folder):
        os.makedirs(folder)


def filesFromFolder(folder, regexp=None):
    contents = os.listdir(folder)
    files = [x for x in contents if regexp.match(x)] if regexp else contents
    return files