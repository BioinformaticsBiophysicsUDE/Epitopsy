# So 12. Nov 19:37:13 CET 2017

import os

def file_ext(filename, lowercase=True, strip_tilde=True):
    '''
    Find the file extension.
    
    :param filename: filename
    :type  filename: str
    :param lowercase: change to lowercase
    :type  lowercase: bool
    :param strip_tilde: remove trailing tilde (hidden file)
    :type  strip_tilde: bool
    :returns: File extension with the leading '.'.
    :returntype: str
    
    Examples::
    
        >>> file_ext('protein_epi.dx')
        '.dx'
        >>> file_ext('protein_epi.PDB~')
        '.pdb'
    
    '''
    ext = os.path.splitext(filename)[1]
    if strip_tilde:
        ext = ext.rstrip('~')
    if lowercase:
        ext = ext.lower()
    return ext


def file_replace_ext(filename, ext):
    '''
    Remove the file extension and add a new one.
    
    :param filename: filename
    :type  filename: str
    :param ext: new extension (with or without leading '.') or empty string
    :type  ext: str or ``None``
    :returns: New filename.
    :returntype: str
    
    Examples::
    
        >>> file_replace_ext('protein_epi.dx', 'pdb')
        'protein_epi.pdb'
        >>> file_replace_ext('protein_epi.dx', None)
        'protein_epi'
        >>> file_replace_ext('protein_epi', 'pdb')
        'protein_epi.pdb'
    
    '''
    if not ext:
        ext = ''
    elif ext[0] != '.':
        ext = '.' + ext
    return os.path.splitext(filename)[0] + ext

