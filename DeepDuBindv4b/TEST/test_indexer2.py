import DuGenCoords

merged_bed_fname = "../ENCODE_DREAM_LINK/train_test/train_regions.blacklistfiltered.merged.bed"
indexer = DuGenCoords.IndexCoords(merged_bed_fname,200,50)

test_index  = 1509126924
print "Test index",test_index
entry = indexer.retrieve_by_index(test_index)
print entry
print entry[0:3]

ret_index = indexer.index_from_coords(entry[0],entry[1],entry[2])
print ret_index

ret_index2 = indexer.index_from_coords('chr22',34260575305, 34260575345)
print ret_index2

print "Identify cell_index",indexer.identify_cell_index(test_index)

tup = indexer.retrieve_by_index(test_index)
tup = '\t'.join([tup[0],str(tup[1]),str(tup[2]),str(tup[3])])
print "Retrieve via index",tup
