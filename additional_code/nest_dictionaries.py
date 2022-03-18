""" This module contains functions that can be used to nest dictionaries. """

def select_nest(dict_2,keys,last,counter=0):
    if last == 0:
        return dict_2
    if counter == 0:
        counter += 1
    if counter == last:
        return dict_2[keys[counter-1]]
    else:
        return select_nest(dict_2[keys[counter-1]],keys,last,counter+1)


def add_dict_branch(dict_1,keys,end_type):
    for count, key in enumerate(keys):
        if count < len(keys)-1:
            if key not in select_nest(dict_1,keys,count).keys():
                select_nest(dict_1,keys,count)[key] = {}
        else:
            if key not in select_nest(dict_1,keys,count).keys():
                select_nest(dict_1,keys,count)[key] = end_type