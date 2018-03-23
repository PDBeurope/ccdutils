def find_element_in_list(array, element):
    """Finds an element in an array. Does not crash if not found

    Arguments:
        array {list} -- array to be searched
        element {any} -- element to be found

    Returns:
        [int?] -- index of element or None
    """
    try:
        index = array.index(element)
        return index
    except ValueError:
        return None
