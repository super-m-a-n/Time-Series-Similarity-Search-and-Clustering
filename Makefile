CXX=g++
CXXFLAGS =-Wall -Wextra -Wno-unused-parameter -Wno-format-truncation -std=c++11
CXXFLAGS += -I. -I./common -I./lsh_folder -I./hypercube_folder #-I./cluster_folder
OBJS_FOLDER = ./objects 
OBJ_COMMON =  ./common/object.o ./common/assist_functions.o ./common/h_hash.o ./common/dataset.o ./common/input_check.o
OBJ_LSH = ./lsh_folder/g_hash.o ./lsh_folder/hash.o ./lsh_folder/lsh_struct.o
OBJ_HYPERCUBE = ./hypercube_folder/f_hash.o ./hypercube_folder/hypercube_class.o
#OBJ_CLUSTER = ./cluster_folder/cluster_info.o 
PROGRAMS = search

all: $(PROGRAMS) mv_objs

target1: search

mv_objs:
	mkdir -p $(OBJS_FOLDER)
	#mv -f $(OBJ_CLUSTER) $(OBJS_FOLDER) 2>/dev/null; true
	mv -f $(OBJ_HYPERCUBE) $(OBJS_FOLDER) 2>/dev/null; true
	mv -f $(OBJ_LSH) $(OBJS_FOLDER) 2>/dev/null; true
	mv -f $(OBJ_COMMON) $(OBJS_FOLDER) 2>/dev/null; true
	mv -f search.o $(OBJS_FOLDER) 2>/dev/null; true

search: $(OBJ_COMMON) $(OBJ_HYPERCUBE) $(OBJ_LSH) search.o 
	$(CXX) $(CXXFLAGS) -o search search.o $(OBJ_COMMON) $(OBJ_HYPERCUBE) $(OBJ_LSH)

.PHONY: clean

clean:
	rm -rf *.o search output* $(OBJS_FOLDER)


###############################################################################################################

cluster_s_c: target3
	./cluster -i ./data/input_small_id -c ./data/cluster.config -o output_s_cl -complete -m Classic

cluster_b_c: target3
	./cluster -i ./data/input_b_id -c ./data/cluster.config -o output_b_cl -complete -m Classic

cluster_s_lsh: target3
	./cluster -i ./data/input_small_id -c ./data/cluster.config -o output_s_cl -complete -m LSH

cluster_b_lsh: target3
	./cluster -i ./data/input_b_id -c ./data/cluster.config -o output_b_cl -complete -m LSH

cluster_s_cube: target3
	./cluster -i ./data/input_small_id -c ./data/cluster.config -o output_s_cl -complete -m Hypercube

cluster_b_cube: target3
	./cluster -i ./data/input_b_id -c ./data/cluster.config -o output_b_cl -complete -m Hypercube

lsh_s_s: target1
	./lsh -i ./data/input_small_id -q ./data/query_small_id -o output_s_s_lsh

lsh_s_b: target1
	./lsh -i ./data/input_small_id -q ./data/query_b_id -o output_s_b_lsh

lsh_b_s: target1
	./lsh -i ./data/input_b_id -q ./data/query_small_id -o output_b_s_lsh

lsh_b_b: target1
	./lsh -i ./data/input_b_id -q ./data/query_b_id -o output_b_b_lsh


cube_s_s: target2
	./cube -i ./data/input_small_id -q ./data/query_small_id -o output_s_s_cube

cube_s_b: target2
	./cube -i ./data/input_small_id -q ./data/query_b_id -o output_s_b_cube

cube_b_s: target2
	./cube -i ./data/input_b_id -q ./data/query_small_id -o output_b_s_cube

cube_b_b: target2
	./cube -i ./data/input_b_id -q ./data/query_b_id -o output_b_b_cube
