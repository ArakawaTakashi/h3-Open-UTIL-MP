import ctypes
import numpy as np

h3oupf = np.ctypeslib.load_library("libh3oup.so","/home/arakawa/work/2024H3OPEN/h3ou2024/h3-Open-UTIL-MP/h3ou/lib")

#=======+=========+=========+=========+=========+=========+=========+=========+

def get_string_buffer(py_str):
    enc_str = py_str.encode("utf-8")
    return ctypes.create_string_buffer(enc_str)

#=======+=========+=========+=========+=========+=========+=========+=========+

def get_strlen_ptr(py_str):
    str_len = len(py_str)
    return ctypes.byref(ctypes.c_int32(str_len))
                       
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_init(comp_name, config_file_name):

    h3oupf.h3oup_init.argtypes = [
         ctypes.POINTER(ctypes.c_char),
         ctypes.POINTER(ctypes.c_int32),
         ctypes.POINTER(ctypes.c_char),
         ctypes.POINTER(ctypes.c_int32)
         ]
    h3oupf.h3oup_init.restype = ctypes.c_void_p

    comp       = ctypes.create_string_buffer(comp_name.encode("utf-8"))
    comp_len   = len(comp_name)
    config     = ctypes.create_string_buffer(config_file_name.encode("utf-8"))
    config_len = len(config_file_name)
    

    h3oupf.h3oup_init(comp, ctypes.byref(ctypes.c_int32(comp_len)), \
                      config, ctypes.byref(ctypes.c_int32(config_len)))
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_my_rank():

    h3oupf.h3oup_get_my_rank.argtypes = [
        ]
    h3oupf.h3oup_get_my_rank.restype = ctypes.c_int32

    return int(h3oupf.h3oup_get_my_rank())

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_my_size():

    h3oupf.h3oup_get_my_rank.argtypes = [
        ]
    h3oupf.h3oup_get_my_rank.restype = ctypes.c_int32

    return int(h3oupf.h3oup_get_my_size())

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_def_grid(grid_index, comp_name, grid_name, nz):
    """ define grid index 
    """
    h3oupf.h3oup_def_grid.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_def_grid.restype = ctypes.c_void_p

    grid_len = len(grid_index)
    glen_ptr = ctypes.byref(ctypes.c_int32(grid_len))

    comp_buffer = get_string_buffer(comp_name)
    cstr_len = get_strlen_ptr(comp_name)

    grid_buffer = get_string_buffer(grid_name)
    gstr_len   = get_strlen_ptr(grid_name)
    
    nz_ptr = ctypes.byref(ctypes.c_int32(nz))

    h3oupf.h3oup_def_grid(grid_index, glen_ptr, comp_buffer, cstr_len, \
                        grid_buffer, gstr_len, nz_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_end_grid_def():
    h3oupf.h3oup_end_grid_def()

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_def_varp(comp_name, data_name, grid_name, num_of_layer):
    """ define recv data
    """
    h3oupf.h3oup_def_varp.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_def_varp.restype = ctypes.c_void_p
    
    comp_buffer = get_string_buffer(comp_name)
    cstr_len    = get_strlen_ptr(comp_name)
    data_buffer = get_string_buffer(data_name)
    dstr_len    = get_strlen_ptr(data_name)
    grid_buffer = get_string_buffer(grid_name)
    gstr_len    = get_strlen_ptr(grid_name)
    
    nz_ptr = ctypes.byref(ctypes.c_int32(num_of_layer))

    h3oupf.h3oup_def_varp(comp_buffer, cstr_len, data_buffer, dstr_len, \
                        grid_buffer, gstr_len, nz_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_def_varg(comp_name, data_name, grid_name, num_of_layer,
                 send_comp_name, send_data_name, recv_mode, intvl, time_lag,
                 mapping_tag, exchange_tag):
    """ define recv data
    """
    h3oupf.h3oup_def_varg.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_def_varp.restype = ctypes.c_void_p
    
    comp_buffer = get_string_buffer(comp_name)
    cstr_len    = get_strlen_ptr(comp_name)
    data_buffer = get_string_buffer(data_name)
    dstr_len    = get_strlen_ptr(data_name)
    grid_buffer = get_string_buffer(grid_name)
    gstr_len    = get_strlen_ptr(grid_name)
    nz_ptr      = ctypes.byref(ctypes.c_int32(num_of_layer))
    sc_buffer   = get_string_buffer(send_comp_name)
    scstr_len   = get_strlen_ptr(send_comp_name)
    sd_buffer   = get_string_buffer(send_data_name)
    sdstr_len   = get_strlen_ptr(send_data_name)
    rm_buffer   = get_string_buffer(recv_mode)
    rmstr_len   = get_strlen_ptr(recv_mode)
    intvl_ptr   = ctypes.byref(ctypes.c_int32(intvl))
    lag_ptr     = ctypes.byref(ctypes.c_int32(time_lag))
    mtag_ptr    = ctypes.byref(ctypes.c_int32(mapping_tag))
    etag_ptr    = ctypes.byref(ctypes.c_int32(exchange_tag))

    h3oupf.h3oup_def_varg(comp_buffer, cstr_len, data_buffer, dstr_len, \
                        grid_buffer, gstr_len, nz_ptr,
                        sc_buffer, scstr_len, sd_buffer, sdstr_len, \
                        rm_buffer, rmstr_len,
                        intvl_ptr, lag_ptr, mtag_ptr, etag_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_end_var_def():
    """ finish variable definition
    """
    h3oupf.h3oup_end_var_def()

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_num_of_put_data():
    h3oupf.h3oup_get_num_of_put_data.argtypes = [
         ctypes.POINTER(ctypes.c_int32)
         ]
    h3oupf.h3oup_get_num_of_put_data.restype = ctypes.c_void_p

    num_of_data = 0
    num_of_data = ctypes.c_int32(num_of_data)

    h3oupf.h3oup_get_num_of_put_data(ctypes.byref(num_of_data))

    return num_of_data.value

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_put_data_name(data_num):
    h3oupf.h3oup_get_put_data_name.argtypes = [
         ctypes.POINTER(ctypes.c_int32),
         ctypes.POINTER(ctypes.c_char),
         ctypes.POINTER(ctypes.c_int32),
    ]
    h3oupf.h3oup_get_put_data_name.restype = ctypes.c_void_p

    dn_ptr    = data_num + 1 # for fortran 
    dn_ptr    = ctypes.c_int32(dn_ptr)
    data_name = "                                        "
    name_buffer = get_string_buffer(data_name)
    nstr_len    = get_strlen_ptr(data_name)
    
    h3oupf.h3oup_get_put_data_name(ctypes.byref(dn_ptr), name_buffer, nstr_len)

    return name_buffer.value.decode().strip(" ")

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_num_of_get_data():
    h3oupf.h3oup_get_num_of_get_data.argtypes = [
         ctypes.POINTER(ctypes.c_int32)
         ]
    h3oupf.h3oup_get_num_of_get_data.restype = ctypes.c_void_p

    num_of_data = 0
    num_of_data = ctypes.c_int32(num_of_data)

    h3oupf.h3oup_get_num_of_get_data(ctypes.byref(num_of_data))

    return num_of_data.value

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_get_data_name(data_num):
    h3oupf.h3oup_get_get_data_name.argtypes = [
         ctypes.POINTER(ctypes.c_int32),
         ctypes.POINTER(ctypes.c_char),
         ctypes.POINTER(ctypes.c_int32),
    ]
    h3oupf.h3oup_get_get_data_name.restype = ctypes.c_void_p

    dn_ptr    = data_num + 1 # for fortran 
    dn_ptr    = ctypes.c_int32(dn_ptr)
    data_name = "                                        "
    name_buffer = get_string_buffer(data_name)
    nstr_len    = get_strlen_ptr(data_name)
    
    h3oupf.h3oup_get_get_data_name(ctypes.byref(dn_ptr), name_buffer, nstr_len)

    return name_buffer.value.decode().strip(" ")

    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_vlayer(data_name):
    h3oupf.h3oup_get_vlayer.argtypes = [
         ctypes.POINTER(ctypes.c_char),
         ctypes.POINTER(ctypes.c_int32),
         ctypes.POINTER(ctypes.c_int32),
    ]
    h3oupf.h3oup_get_vlayer.restype = ctypes.c_void_p

    data_buffer = get_string_buffer(data_name)
    dstr_len    = get_strlen_ptr(data_name)
    vlayer    = 1
    vlayer    = ctypes.c_int32(vlayer)
    
    h3oupf.h3oup_get_vlayer(data_buffer, dstr_len, ctypes.byref(vlayer))

    return vlayer.value

    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_set_interpolation_table(my_name, send_comp, send_grid, recv_comp, recv_grid,
                                mapping_tag, send_index, recv_index, coef):
    """ define mapping table index
    """
    h3oupf.h3oup_set_interpolation_table.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # my_name
        ctypes.POINTER(ctypes.c_int32),         # my_name length
        ctypes.POINTER(ctypes.c_char),          # send_comp
        ctypes.POINTER(ctypes.c_int32),         # send_comp length
        ctypes.POINTER(ctypes.c_char),          # send_grid
        ctypes.POINTER(ctypes.c_int32),         # send_grid length
        ctypes.POINTER(ctypes.c_char),          # recv_comp
        ctypes.POINTER(ctypes.c_int32),         # recv_comp length
        ctypes.POINTER(ctypes.c_char),          # recv_grid
        ctypes.POINTER(ctypes.c_int32),         # recv_grid length
        ctypes.POINTER(ctypes.c_int32),         # mapping_tag
        np.ctypeslib.ndpointer(dtype=np.int32), # send_index
        np.ctypeslib.ndpointer(dtype=np.int32), # recv_index
        np.ctypeslib.ndpointer(dtype=np.float64),# coefficient
        ctypes.POINTER(ctypes.c_int32)          # size of index array
        ]
    h3oupf.h3oup_set_interpolation_table.restype = ctypes.c_void_p

    myname_buffer = get_string_buffer(my_name)
    mnstr_len     = get_strlen_ptr(my_name)
    scomp_buffer  = get_string_buffer(send_comp)
    scstr_len     = get_strlen_ptr(send_comp)
    sgrid_buffer  = get_string_buffer(send_grid)
    sgstr_len     = get_strlen_ptr(send_grid)
    rcomp_buffer  = get_string_buffer(recv_comp)
    rcstr_len     = get_strlen_ptr(recv_comp)
    rgrid_buffer  = get_string_buffer(recv_grid)
    rgstr_len     = get_strlen_ptr(recv_grid)

    mtag_ptr      = ctypes.byref(ctypes.c_int32(mapping_tag))
    grid_len      = len(send_index)
    glen_ptr      = ctypes.byref(ctypes.c_int32(grid_len))
    
    h3oupf.h3oup_set_interpolation_table(myname_buffer, mnstr_len,
                                 scomp_buffer, scstr_len, sgrid_buffer, sgstr_len,
                                 rcomp_buffer, rcstr_len, rgrid_buffer, rgstr_len,
                                 mtag_ptr, send_index, recv_index, coef, glen_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_set_interpolation_table_no_index(my_name, send_comp, send_grid, \
                                         recv_comp, recv_grid,
                                   mapping_tag):
    """ define mapping table index
    """
    h3oupf.h3oup_set_interpolation_table_no_index.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # my_name
        ctypes.POINTER(ctypes.c_int32),         # my_name length
        ctypes.POINTER(ctypes.c_char),          # send_comp
        ctypes.POINTER(ctypes.c_int32),         # send_comp length
        ctypes.POINTER(ctypes.c_char),          # send_grid
        ctypes.POINTER(ctypes.c_int32),         # send_grid length
        ctypes.POINTER(ctypes.c_char),          # recv_comp
        ctypes.POINTER(ctypes.c_int32),         # recv_comp length
        ctypes.POINTER(ctypes.c_char),          # recv_grid
        ctypes.POINTER(ctypes.c_int32),         # recv_grid length
        ctypes.POINTER(ctypes.c_int32)          # mapping_tag
        ]
    h3oupf.h3oup_set_interpolation_table_no_index.restype = ctypes.c_void_p

    myname_buffer = get_string_buffer(my_name)
    mnstr_len     = get_strlen_ptr(my_name)
    scomp_buffer  = get_string_buffer(send_comp)
    scstr_len     = get_strlen_ptr(send_comp)
    sgrid_buffer  = get_string_buffer(send_grid)
    sgstr_len     = get_strlen_ptr(send_grid)
    rcomp_buffer  = get_string_buffer(recv_comp)
    rcstr_len     = get_strlen_ptr(recv_comp)
    rgrid_buffer  = get_string_buffer(recv_grid)
    rgstr_len     = get_strlen_ptr(recv_grid)
    mtag_ptr      = ctypes.byref(ctypes.c_int32(mapping_tag))

    h3oupf.h3oup_set_interpolation_table_no_index(myname_buffer, mnstr_len,
                                 scomp_buffer, scstr_len, sgrid_buffer, sgstr_len,
                                 rcomp_buffer, rcstr_len, rgrid_buffer, rgstr_len,
                                 mtag_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_set_mapping_table(my_name, send_comp, send_grid, recv_comp, recv_grid,
                          mapping_tag, send_index, recv_index):
    """ define mapping table index
    """
    h3oupf.h3oup_set_mapping_table.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # my_name
        ctypes.POINTER(ctypes.c_int32),         # my_name length
        ctypes.POINTER(ctypes.c_char),          # send_comp
        ctypes.POINTER(ctypes.c_int32),         # send_comp length
        ctypes.POINTER(ctypes.c_char),          # send_grid
        ctypes.POINTER(ctypes.c_int32),         # send_grid length
        ctypes.POINTER(ctypes.c_char),          # recv_comp
        ctypes.POINTER(ctypes.c_int32),         # recv_comp length
        ctypes.POINTER(ctypes.c_char),          # recv_grid
        ctypes.POINTER(ctypes.c_int32),         # recv_grid length
        ctypes.POINTER(ctypes.c_int32),         # mapping_tag
        np.ctypeslib.ndpointer(dtype=np.int32), # send_index
        np.ctypeslib.ndpointer(dtype=np.int32), # recv_index
        ctypes.POINTER(ctypes.c_int32)          # size of index array
        ]
    h3oupf.h3oup_set_mapping_table.restype = ctypes.c_void_p

    myname_buffer = get_string_buffer(my_name)
    mnstr_len     = get_strlen_ptr(my_name)
    scomp_buffer  = get_string_buffer(send_comp)
    scstr_len     = get_strlen_ptr(send_comp)
    sgrid_buffer  = get_string_buffer(send_grid)
    sgstr_len     = get_strlen_ptr(send_grid)
    rcomp_buffer  = get_string_buffer(recv_comp)
    rcstr_len     = get_strlen_ptr(recv_comp)
    rgrid_buffer  = get_string_buffer(recv_grid)
    rgstr_len     = get_strlen_ptr(recv_grid)

    mtag_ptr      = ctypes.byref(ctypes.c_int32(mapping_tag))
    grid_len      = len(send_index)
    glen_ptr      = ctypes.byref(ctypes.c_int32(grid_len))
    
    h3oupf.h3oup_set_mapping_table(myname_buffer, mnstr_len,
                                 scomp_buffer, scstr_len, sgrid_buffer, sgstr_len,
                                 rcomp_buffer, rcstr_len, rgrid_buffer, rgstr_len,
                                 mtag_ptr, send_index, recv_index, glen_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_set_mapping_table_no_index(my_name, send_comp, send_grid, \
                                   recv_comp, recv_grid,
                                   mapping_tag):
    """ define mapping table index
    """
    h3oupf.h3oup_set_mapping_table_no_index.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # my_name
        ctypes.POINTER(ctypes.c_int32),         # my_name length
        ctypes.POINTER(ctypes.c_char),          # send_comp
        ctypes.POINTER(ctypes.c_int32),         # send_comp length
        ctypes.POINTER(ctypes.c_char),          # send_grid
        ctypes.POINTER(ctypes.c_int32),         # send_grid length
        ctypes.POINTER(ctypes.c_char),          # recv_comp
        ctypes.POINTER(ctypes.c_int32),         # recv_comp length
        ctypes.POINTER(ctypes.c_char),          # recv_grid
        ctypes.POINTER(ctypes.c_int32),         # recv_grid length
        ctypes.POINTER(ctypes.c_int32),         # mapping_tag
        ]
    h3oupf.h3oup_set_mapping_table_no_index.restype = ctypes.c_void_p

    myname_buffer = get_string_buffer(my_name)
    mnstr_len     = get_strlen_ptr(my_name)
    scomp_buffer  = get_string_buffer(send_comp)
    scstr_len     = get_strlen_ptr(send_comp)
    sgrid_buffer  = get_string_buffer(send_grid)
    sgstr_len     = get_strlen_ptr(send_grid)
    rcomp_buffer  = get_string_buffer(recv_comp)
    rcstr_len     = get_strlen_ptr(recv_comp)
    rgrid_buffer  = get_string_buffer(recv_grid)
    rgstr_len     = get_strlen_ptr(recv_grid)
    mtag_ptr      = ctypes.byref(ctypes.c_int32(mapping_tag))

    h3oupf.h3oup_set_mapping_table_no_index(myname_buffer, mnstr_len,
                                 scomp_buffer, scstr_len, sgrid_buffer, sgstr_len,
                                 rcomp_buffer, rcstr_len, rgrid_buffer, rgstr_len,
                                 mtag_ptr)
    
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_init_time(time_array):
    """ set initial integration time
    
    """
    h3oupf.h3oup_init_time.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), # time array
        ctypes.POINTER(ctypes.c_int32)          # size of time array
        ]
    h3oupf.h3oup_init_time.restype = ctypes.c_void_p
    
    array_len      = len(time_array)
    glen_ptr      = ctypes.byref(ctypes.c_int32(array_len))
    
    h3oupf.h3oup_init_time(time_array, glen_ptr)

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_set_time(my_name, time_array, delta_t):
    """ set initial integration time
    
    """
    h3oupf.h3oup_set_time.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # my_name
        ctypes.POINTER(ctypes.c_int32),         # my_name length
        np.ctypeslib.ndpointer(dtype=np.int32), # time array
        ctypes.POINTER(ctypes.c_int32),         # size of time array
        ctypes.POINTER(ctypes.c_int32)          # delta_t
        ]
    h3oupf.h3oup_set_time.restype = ctypes.c_void_p

    myname_buffer = get_string_buffer(my_name)
    mnstr_len     = get_strlen_ptr(my_name)
    array_len      = len(time_array)
    glen_ptr      = ctypes.byref(ctypes.c_int32(array_len))
    delt_ptr      = ctypes.byref(ctypes.c_int32(delta_t))
    
    h3oupf.h3oup_set_time(myname_buffer, mnstr_len, time_array, glen_ptr, delt_ptr)

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_put_data_1d(var_name, data):
    """ put 1d data
    
    """
    h3oupf.h3oup_put_data_1d.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # data_name
        ctypes.POINTER(ctypes.c_int32),         # data_name length
        np.ctypeslib.ndpointer(dtype=np.float64), # data array
        ctypes.POINTER(ctypes.c_int32)            # size of data array
        ]
    h3oupf.h3oup_put_data_1d.restype = ctypes.c_void_p

    varname_buffer = get_string_buffer(var_name)
    varstr_len     = get_strlen_ptr(var_name)
    array_len = len(data)
    dlen_ptr  = ctypes.byref(ctypes.c_int32(array_len))

    h3oupf.h3oup_put_data_1d(varname_buffer, varstr_len, data, dlen_ptr)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_put_data_25d(var_name, data):
    """ put 25d data
    
    """
    h3oupf.h3oup_put_data_25d.argtypes = [
        ctypes.POINTER(ctypes.c_char),          # data_name
        ctypes.POINTER(ctypes.c_int32),         # data_name length
        np.ctypeslib.ndpointer(dtype=np.float64), # data array
        ctypes.POINTER(ctypes.c_int32),           # size of data array
        ctypes.POINTER(ctypes.c_int32)            # size of data array
        ]
    h3oupf.h3oup_put_data_25d.restype = ctypes.c_void_p

    varname_buffer = get_string_buffer(var_name)
    varstr_len     = get_strlen_ptr(var_name)
    array_len1 = data.shape[0]
    dlen_ptr1  = ctypes.byref(ctypes.c_int32(array_len1))
    array_len2 = data.shape[1]
    dlen_ptr2  = ctypes.byref(ctypes.c_int32(array_len2))

    h3oupf.h3oup_put_data_25d(varname_buffer, varstr_len, data, dlen_ptr1, dlen_ptr2)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_data_1d(var_name, data):
    """ get 1d data
    
    """
    h3oupf.h3oup_get_data_1d.argtypes = [
        ctypes.POINTER(ctypes.c_char),            # data_name
        ctypes.POINTER(ctypes.c_int32),           # data_name length
        np.ctypeslib.ndpointer(dtype=np.float64), # data array
        ctypes.POINTER(ctypes.c_int32),           # size of data array
        ctypes.POINTER(ctypes.c_int32)            # recv flag
        ]

    h3oupf.h3oup_get_data_1d.restype = ctypes.c_void_p

    varname_buffer = get_string_buffer(var_name)
    varstr_len     = get_strlen_ptr(var_name)
    array_len = len(data)
    dlen_ptr  = ctypes.byref(ctypes.c_int32(array_len))

    is_get_ok = True
    ok_flag    = 1
    
    h3oupf.h3oup_get_data_1d(varname_buffer, varstr_len, data, dlen_ptr, \
                             ctypes.byref(ctypes.c_int32(ok_flag)))

    if (ok_flag == 0):
        is_get_ok = False

    return is_get_ok

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_get_data_25d(var_name, data):
    """ get 25d data
    
    """
    h3oupf.h3oup_get_data_25d.argtypes = [
        ctypes.POINTER(ctypes.c_char),            # data_name
        ctypes.POINTER(ctypes.c_int32),           # data_name length
        np.ctypeslib.ndpointer(dtype=np.float64), # data array
        ctypes.POINTER(ctypes.c_int32),           # size of data array
        ctypes.POINTER(ctypes.c_int32),           # size of data array
        ctypes.POINTER(ctypes.c_int32)            # recv flag
        ]

    h3oupf.h3oup_get_data_25d.restype = ctypes.c_void_p

    varname_buffer = get_string_buffer(var_name)
    varstr_len     = get_strlen_ptr(var_name)
    array_len1 = data.shape[0]
    dlen_ptr1  = ctypes.byref(ctypes.c_int32(array_len1))
    array_len2 = data.shape[1]
    dlen_ptr2  = ctypes.byref(ctypes.c_int32(array_len2))

    is_get_ok = True
    ok_flag    = 1
    
    h3oupf.h3oup_get_data_25d(varname_buffer, varstr_len, data, dlen_ptr1, dlen_ptr2, \
                              ctypes.byref(ctypes.c_int32(ok_flag)))

    if (ok_flag == 0):
        is_get_ok = False

    return is_get_ok
        
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_coupling_end(time_array):
    """ finalize coupling
    
    """
    h3oupf.h3oup_coupling_end.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.int32), # time array
        ctypes.POINTER(ctypes.c_int32),         # size of time array
        ]
    h3oupf.h3oup_coupling_end.restype = ctypes.c_void_p

    array_len      = len(time_array)
    glen_ptr      = ctypes.byref(ctypes.c_int32(array_len))

    h3oupf.h3oup_coupling_end(time_array, glen_ptr)

#=======+=========+=========+=========+=========+=========+=========+=========+
#=======+=========+=========+=========+=========+=========+=========+=========+
#=======+=========+=========+=========+=========+=========+=========+=========+

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_end():
    """ finalize coupling
    
    """
    h3oupf.h3oup_end.restype = ctypes.c_void_p

    h3oupf.h3oup_end()

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_scalar(target_name, val):

    h3oupf.h3oup_send_int_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_int_scalar.restype = ctypes.c_void_p

    h3oupf.h3oup_send_real_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float)
        ]
    h3oupf.h3oup_send_real_scalar.restype = ctypes.c_void_p
    
    h3oupf.h3oup_send_double_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double)
        ]
    h3oupf.h3oup_send_double_scalar.restype = ctypes.c_void_p
    
    comp_str      = ctypes.create_string_buffer(target_name.encode("utf-8"))
    comp_len      = len(target_name)

    if isinstance(val, int):
        h3oupf.h3oup_send_int_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(val)))
    elif isinstance(val, float):
        h3oupf.h3oup_send_real_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_float(val)))
    else:
        h3oupf.h3oup_send_double_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_double(val)))
            
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_scalar(target_name, val):

    h3oupf.h3oup_recv_int_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_int_scalar.restype = ctypes.c_void_p

    h3oupf.h3oup_recv_real_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float)
        ]
    h3oupf.h3oup_recv_real_scalar.restype = ctypes.c_void_p
    
    h3oupf.h3oup_recv_double_scalar.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double)
        ]
    h3oupf.h3oup_recv_double_scalar.restype = ctypes.c_void_p
    
    comp_str      = ctypes.create_string_buffer(target_name.encode("utf-8"))
    comp_len      = len(target_name)

    if isinstance(val, int):
        int_val = ctypes.c_int32(val)
        h3oupf.h3oup_recv_int_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(int_val))
        return int_val.value
    elif isinstance(val, float):
        float_val = ctypes.c_float(val)
        h3oupf.h3oup_recv_real_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(float_val))
        return float_val.value 
    else:
        double_val = ctypes.c_doube(val)
        h3oupf.h3oup_recv_double_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(double_val))
        return double_val.value 
            
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_array(target_name, val):

    h3oupf.h3oup_send_int_array.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_int_array.restype = ctypes.c_void_p

    h3oupf.h3oup_send_real_array.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_real_array.restype = ctypes.c_void_p
    
    h3oupf.h3oup_send_double_array.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_double_array.restype = ctypes.c_void_p
    
    comp_str      = ctypes.create_string_buffer(target_name.encode("utf-8"))
    comp_len      = len(target_name)
    array_len     = len(val)
    
    if isinstance(val[0], int):
        c_array = (ctypes.c_int32 * len(val))(*val)
        h3oupf.h3oup_send_int_array(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                    ctypes.byref(ctypes.c_int32(array_len)))
    elif isinstance(va[0], float):
        c_array = (ctypes.c_float * len(val))(*val)
        h3oupf.h3oup_send_real_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.cast(c_array, ctype.POINTER(ctypes.c_float)), \
                                      ctypes.byref(ctypes.c_int32(array_len)))
    else:
        c_array = (ctypes.c_double * len(val))(*val)
        h3oupf.h3oup_send_real_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.cast(c_array, ctype.POINTER(ctypes.c_double)), \
                                      ctypes.byref(ctypes.c_int32(array_len)))
            
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_array(target_name, val):

    h3oupf.h3oup_recv_int_array.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_int_array.restype = ctypes.c_void_p

    h3oupf.h3oup_recv_real_array.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_real_array.restype = ctypes.c_void_p
    
    h3oupf.h3oup_recv_double_array.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_double_array.restype = ctypes.c_void_p
    
    comp_str      = ctypes.create_string_buffer(target_name.encode("utf-8"))
    comp_len      = len(target_name)
    array_len     = len(val)
    
    if isinstance(val[0], int):
        c_array = (ctypes.c_int32 * len(val))(*val)
        h3oupf.h3oup_recv_int_array(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                    ctypes.byref(ctypes.c_int32(array_len)))
        return list(c_array)
    
    elif isinstance(va[0], float):
        c_array = (ctypes.c_float * len(val))(*val)
        h3oupf.h3oup_recv_real_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.cast(c_array, ctype.POINTER(ctypes.c_float)), \
                                      ctypes.byref(ctypes.c_int32(array_len)))
        return list(c_array)

    else:
        c_array = (ctypes.c_double * len(val))(*val)
        h3oupf.h3oup_recv_real_scalar(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.cast(c_array, ctype.POINTER(ctypes.c_double)), \
                                      ctypes.byref(ctypes.c_int32(array_len)))
        return list(c_array)
            
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_model_int(target_name, target_pe, data):

    h3oupf.h3oup_send_model_int.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_model_int.restype = ctypes.c_void_p
    
    comp_str = ctypes.create_string_buffer(target_name.encode("utf-8"))
    comp_len = len(target_name)
    data_len = len(data)

    c_array = (ctypes.c_int32 * len(data))(*data)
    h3oupf.h3oup_send_model_int(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(target_pe)), \
                                ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                ctypes.byref(ctypes.c_int32(data_len)))

    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_model_real(target_name, target_pe, data):

    h3oupf.h3oup_send_model_real.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_model_real.restype = ctypes.c_void_p
    
    c_array = (ctypes.c_float * len(data))(*data)
    h3oupf.h3oup_send_model_real(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(target_pe)), \
                                 ctypes.cast(c_array, ctypes.POINTER(ctypes.c_float)), \
                                 ctypes.byref(ctypes.c_int32(data_len)))
     
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_model_double(target_name, target_pe, data):

    h3oupf.h3oup_send_model_double.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_model_double.restype = ctypes.c_void_p

    comp_str = ctypes.create_string_buffer(target_name.encode("utf-8"))
    comp_len = len(target_name)
    data_len = len(data)

    c_array = (ctypes.c_double * len(data))(*data)
    h3oupf.h3oup_send_model_double(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(target_pe)), \
                                   ctypes.cast(c_array, ctypes.POINTER(ctypes.c_double)), \
                                   ctypes.byref(ctypes.c_int32(data_len)))
    
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_model_int(source_name, source_pe, data):

    h3oupf.h3oup_recv_model_int.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_model_int.restype = ctypes.c_void_p
    
    comp_str = ctypes.create_string_buffer(source_name.encode("utf-8"))
    comp_len = len(source_name)
    data_len = len(data)

    c_array = (ctypes.c_int32 * len(data))(*data)
    h3oupf.h3oup_recv_model_int(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(source_pe)), \
                                ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_model_real(source_name, source_pe, data):

    h3oupf.h3oup_recv_model_real.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_model_real.restype = ctypes.c_void_p
    
    comp_str = ctypes.create_string_buffer(source_name.encode("utf-8"))
    comp_len = len(source_name)
    data_len = len(data)

    c_array = (ctypes.c_float * len(data))(*data)
    h3oupf.h3oup_recv_model_real(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(source_pe)), \
                                 ctypes.cast(c_array, ctypes.POINTER(ctypes.c_float)), \
                                 ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)

    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_model_double(source_name, source_pe, data):

    h3oupf.h3oup_recv_model_double.argtypes = [
        ctypes.POINTER(ctypes.c_char),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_model_double.restype = ctypes.c_void_p
    
    comp_str = ctypes.create_string_buffer(source_name.encode("utf-8"))
    comp_len = len(source_name)
    data_len = len(data)

    c_array = (ctypes.c_double * len(data))(*data)
    h3oupf.h3oup_recv_model_double(comp_str, ctypes.byref(ctypes.c_int32(comp_len)), ctypes.byref(ctypes.c_int32(source_pe)), \
                                   ctypes.cast(c_array, ctypes.POINTER(ctypes.c_double)), \
                                   ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_bcast_local_int(source_pe, data):
    h3oupf.h3oup_bcast_local_int.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_bcast_local_int.restype = ctypes.c_void_p

    data_len = len(data)

    c_array = (ctypes.c_int * len(data))(*data)

    h3oupf.h3oup_bcast_local_int(ctypes.byref(ctypes.c_int32(source_pe)), \
                                 ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                 ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_bcast_local_real(source_pe, data):
    h3oupf.h3oup_bcast_local_int.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_bcast_local_real.restype = ctypes.c_void_p

    data_len = len(data)

    c_array = (ctypes.c_float * len(data))(*data)

    h3oupf.h3oup_bcast_local_real(ctypes.byref(ctypes.c_int32(source_pe)), \
                                  ctypes.cast(c_array, ctypes.POINTER(ctypes.c_float)), \
                                  ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_bcast_local_double(source_pe, data):
    h3oupf.h3oup_bcast_local_int.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_bcast_local_double.restype = ctypes.c_void_p

    data_len = len(data)

    c_array = (ctypes.c_double * len(data))(*data)

    h3oupf.h3oup_bcast_local_int(ctypes.byref(ctypes.c_int32(source_pe)), \
                                 ctypes.cast(c_array, ctypes.POINTER(ctypes.c_double)), \
                                 ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_local_int(target_pe, data):

    h3oupf.h3oup_send_local_int.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_local_int.restype = ctypes.c_void_p
    
    data_len = len(data)

    c_array = (ctypes.c_int32 * len(data))(*data)
    h3oupf.h3oup_send_local_int(ctypes.byref(ctypes.c_int32(target_pe)), \
                                ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                ctypes.byref(ctypes.c_int32(data_len)))

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_local_real(target_pe, data):

    h3oupf.h3oup_send_local_real.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_local_real.restype = ctypes.c_void_p
    
    data_len = len(data)

    c_array = (ctypes.c_float * len(data))(*data)
    h3oupf.h3oup_send_local_real(ctypes.byref(ctypes.c_int32(target_pe)), \
                                 ctypes.cast(c_array, ctypes.POINTER(ctypes.c_float)), \
                                 ctypes.byref(ctypes.c_int32(data_len)))

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_send_local_double(target_pe, data):

    h3oupf.h3oup_send_local_double.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_send_local_double.restype = ctypes.c_void_p
    
    data_len = len(data)

    c_array = (ctypes.c_double * len(data))(*data)
    h3oupf.h3oup_send_local_double(ctypes.byref(ctypes.c_int32(target_pe)), \
                                   ctypes.cast(c_array, ctypes.POINTER(ctypes.c_double)), \
                                   ctypes.byref(ctypes.c_int32(data_len)))

#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_local_int(source_pe, data):

    h3oupf.h3oup_recv_local_int.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_local_int.restype = ctypes.c_void_p
    
    data_len = len(data)

    c_array = (ctypes.c_int32 * len(data))(*data)
    h3oupf.h3oup_recv_local_int(ctypes.byref(ctypes.c_int32(source_pe)), \
                                ctypes.cast(c_array, ctypes.POINTER(ctypes.c_int)), \
                                ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_local_real(source_pe, data):

    h3oupf.h3oup_recv_local_real.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_float),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_local_real.restype = ctypes.c_void_p
    
    data_len = len(data)

    c_array = (ctypes.c_float * len(data))(*data)
    h3oupf.h3oup_recv_local_real(ctypes.byref(ctypes.c_int32(source_pe)), \
                                 ctypes.cast(c_array, ctypes.POINTER(ctypes.c_float)), \
                                 ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+

def h3ou_recv_local_double(source_pe, data):

    h3oupf.h3oup_recv_local_real.argtypes = [
        ctypes.POINTER(ctypes.c_int32),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int32)
        ]
    h3oupf.h3oup_recv_local_double.restype = ctypes.c_void_p
    
    data_len = len(data)

    c_array = (ctypes.c_double * len(data))(*data)
    h3oupf.h3oup_recv_local_double(ctypes.byref(ctypes.c_int32(source_pe)), \
                                   ctypes.cast(c_array, ctypes.POINTER(ctypes.c_double)), \
                                   ctypes.byref(ctypes.c_int32(data_len)))
    return list(c_array)
    
#=======+=========+=========+=========+=========+=========+=========+=========+
