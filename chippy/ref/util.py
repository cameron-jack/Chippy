
# we ignore mitochondria
mouse_chroms = map(str, range(1,20)+['X', 'Y'])
human_chroms = map(str, range(1,23)+['X', 'Y'])
chroms = dict(mouse=mouse_chroms, human=human_chroms)
