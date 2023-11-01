# import pysam

def fetch_upto_next_ID(bam_reader,r): # r is the current record
    records_list=[r] #the records with the same readID
#     print(records_list[0].to_string())
    id_wanted=r.query_name
    try:
        record_current=bam_reader.__next__()
        id_current=record_current.query_name
        while id_current==id_wanted:

            records_list+=[record_current]
            record_current=bam_reader.__next__()
            id_current=record_current.query_name
        return records_list, record_current
    except StopIteration:
        return records_list, None

def fetch_upto_next_ID_aligned(bam_reader,r): # r is the current record
    records_list=[r] #the records with the same readID
#     print(records_list[0].to_string())
    id_wanted=r.query_name
    try:
        record_current=nextAligned(bam_reader)
        id_current=record_current.query_name
        while id_current==id_wanted:

            records_list+=[record_current]
            record_current=nextAligned(bam_reader)
            id_current=record_current.query_name
        return records_list, record_current
    except StopIteration:
        return records_list, None

def nextAligned(bam_reader):
    goon=True
    while goon:
        r=bam_reader.__next__()
        if not r.is_unmapped:
            goon=False
    return r

def fetch_ID(bam_reader,r,id_wanted):
    while ((not r is None) and r.query_name<id_wanted):
        try:
            r=bam_reader.__next__()
        except StopIteration:
            return [], None
        
    if not r is None and r.query_name==id_wanted:
        return fetch_upto_next_ID(bam_reader,r)
    else:
        return [], r


def _lt_natural(x,y, naturalSort):
    if naturalSort:
        return ([ int(c) if c.isdigit() else c.lower() for c in x.split(':') ] < [ int(c) if c.isdigit() else c.lower() for c in y.split(':') ] )
    else:
        return x<y



