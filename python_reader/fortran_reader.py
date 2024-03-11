import numpy as np 
import struct 

def read_record(f, dtype='f4'):
    # Following the description of the standard, found in:
    #  https://stackoverflow.com/questions/15608421/inconsistent-record-marker-while-reading-fortran-unformatted-file
    # This method works always, even when there are split records (because they are longer than 2^31 bytes)!!!!
    head=struct.unpack('i',f.read(4))
    recl=head[0]

    if recl>0:
        numval=recl//np.dtype(dtype).itemsize
        data=np.fromfile(f,dtype=dtype,count=numval)
        endrec=struct.unpack('i',f.read(4))[0]
        #assert recl==endrec
    else:
        the_bytes=f.read(abs(recl))
        endrec=struct.unpack('i',f.read(4))[0]
        while recl<0:
            head=struct.unpack('i',f.read(4))
            recl=head[0]
            the_bytes=the_bytes+f.read(abs(recl))
            endrec=struct.unpack('i',f.read(4))[0]

        if dtype=='f4':
            dtype2='f'
            len_data=4
        elif dtype=='i4':
            dtype2='i'
            len_data=4
        elif dtype=='f8':
            dtype2='d'
            len_data=4
        elif dtype=='i8':
            dtype2='q'
            len_data=4
        else:
            print('Unknown data type!!!!')
            raise ValueError
        dtype2='{:}'.format(len(the_bytes)//len_data)+dtype2
        #print(dtype2)

        data=struct.unpack(dtype2,the_bytes)

    if len(data)==0:
        return np.array([], dtype=dtype)

    return data

def skip_record(f):
    head=struct.unpack('i',f.read(4))
    recl=head[0]

    if recl>0:
        f.seek(recl,1)
        endrec=struct.unpack('i',f.read(4))[0]
        #assert recl==endrec
    else:
        f.seek(abs(recl),1)
        endrec=struct.unpack('i',f.read(4))[0]
        while recl<0:
            head=struct.unpack('i',f.read(4))
            recl=head[0]
            f.seek(abs(recl),1)
            endrec=struct.unpack('i',f.read(4))[0]

    return