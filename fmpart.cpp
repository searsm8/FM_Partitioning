//fmpart.cpp
//Mark Sears and Navya Sreeram
//Fidducia-Mattheyses method to implement bi-partitioning

#include <stdio.h>
#include <stdlib.h>  
#include <stdint.h>
#include <string.h>
#include <string>
#include <list>
#include <map> 
#include <vector> 
#include <iterator>
#include <iostream>
#include <ctime>

using namespace std;

class NET;
class CELL;

//***************************//
//  DECLARE GLOBAL VARIABLES //
//***************************//

char file_path[] = "./FPGA Benchmarks/FPGA-example"; //path where design files can be found

multimap <int, CELL*> gain_bucket [2];

vector <CELL*> cells; //pointers to all cells in the design

map <string, CELL*> cells_map; //map of cells. Key = cell_name, value = cell pointer
				//useful for accessing cells by name

int cut_size = 0; //total cuts across partitions

float MIN_CUT_RATIO = .4; //partitioning constraint. The minimum amount of area allowed in a partition
			//e.g. 0.4 means at least 40% of total area must be in each partition		
			
int partition_area[2] = {0, 0}; //used to track area of partitions

vector <NET*> nets; //pointers to all nets in the design


int min_cut_size;
int initial_partition_cut_size = 0;
int base_partition = 0;
int pass_count = 1;

bool positive_gains_exhausted = false; //true after the initial gain bucket has been exhausted of positive gains
						//used for updating the best_partition
std::clock_t start_time;
std::clock_t timestamp;


//***************************//
//     FUNCTION HEADERS	     //
//***************************//
void readCellsFile();
void readNetsFile();
void computeGains();
vector <CELL*> getNeighbors(CELL* current_cell);
void initialPartition();
void printPartitionAreas();
void findCutsize();
bool FMPartitionPass();
bool swapCells(CELL* cell1, CELL* cell2);
void printTimestamp();
bool performNextMove();
bool cellCanMove(CELL* base_cell);
int printCutsize();

class NET {
	public:
	char name[100];
	int num_pins; //total number of pins on this net
	bool cutstate; //true if the net is cut across the partition
	bool critical; //a net if critical if there is a cell move that affects cutstate
	bool has_locked[2];//a net is "dead" if no legal cell swaps remain that could change the cutstate
			//this occurs when there is a locked or fixed cell in both partitions on this net
			
	int partition_count[2]; //count of how many cells in each partition on this net
	
	vector <CELL*> cell_list; //pointers to all cells on this net			
};

class CELL {
	public:
	char name[100];
	char type[20];
	int area; 	//physical area of the cell
	bool fixed; 	//true for cells that can never be moved
	bool locked; 	//true for cells locked for the rest of the pass
	int gain; 	//gain is the change in cutset if this cell were to swap partitions
	int partition; 	//which partition this cell is in: 0 or 1
	int best_partition; //used to remember this cell's partition in the best solution found so far
	multimap <int, CELL*> :: iterator gain_itr; //an iterator that points to this cell's location in the gain buckets		
	
	vector <NET*> net_list; //pointers to all nets this cell is on
		
	
	//change the partition of the cell. Update gains for this cell and all neighbors. Update cutsize.
	void changePartition()
	{		
		//remove the cell from the bucket i.e. "lock" the cell until the next pass
		gain_bucket[partition].erase(gain_itr);				
		locked = true;	
		
		//update cut_size. cutsize is reduced by the gain of the moved cell.
		cut_size -= gain;		
			
		int base_partition = partition;
		
		//change partition for this cell
		partition_area[partition] -= area; //decrease area in old partition
		partition = !partition;		
		partition_area[partition] += area; //increase area in new partition
		
		if(gain < 0)
			positive_gains_exhausted = true;
											
		//Update the gain of all neighboring cells		
		//Loop through all nets the cell is on. Update gains for neighboring cells on critical nets
		for(int i = 0; i < net_list.size(); i++)
		{						
			//if the net has a locked (or fixed) cell in both partitions,
			//then the net is "dead" and cannot cause any change in cutstate (prevents wasting time)
			if(net_list[i]->has_locked[0] && net_list[i]->has_locked[1])
			{				
				continue;
			}
							
			//create "From" and "To" partition counts
			int from = net_list[i]->partition_count[!partition];   //i.e. F(n), where the cell came from
			int to   = net_list[i]->partition_count[partition];    //i.e. T(n), where the cell is going to
			
		//check critical nets before the move
			if(to == 0)
			{
			//increment gain on all free cells on this net
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if((!net_list[i]->cell_list[k]->locked))		
					{
						net_list[i]->cell_list[k]->updateGain(1);
					}
				}
			}								
			else if (to == 1)
			{
			//decrement gain of the only "to" cell
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if(!net_list[i]->cell_list[k]->locked && base_partition != net_list[i]->cell_list[k]->partition)
					{						
						net_list[i]->cell_list[k]->updateGain(-1);
						break;
					}
				}
			}
							
			//update from and to counts
			from--;
			net_list[i]->partition_count[!partition]--;
			to++;
			net_list[i]->partition_count[partition]++;
			
			//check critical nets after the move
			if(from == 0)
			{ 
			//decrement gain on all free cells on this net
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if(!net_list[i]->cell_list[k]->locked)		
					{
						net_list[i]->cell_list[k]->updateGain(-1);
					}
				}
			}								
			else if (from == 1)
			{
			//increment gain of the only "from" cell
				for(int k = 0; k < net_list[i]->cell_list.size(); k++)
				{
					if(!net_list[i]->cell_list[k]->locked && base_partition == net_list[i]->cell_list[k]->partition)
					{
						net_list[i]->cell_list[k]->updateGain(1);
						break;
					}
				}
			}
			
			//update the has_locked member for this net for the "To" partition
			net_list[i]->has_locked[!base_partition] = true;
		}	
	}
	
	
	void updateGain(int gain_change)
	{		
		//if the cell is locked, do nothing. Only update free cells.
		//because locked cells are no longer considered for this pass
		if(locked || fixed) return;
		
		//remove the old entry in gain bucket
		gain_bucket[partition].erase(gain_itr);
		
		gain += gain_change;
						
		//add the new entry to gain bucket, and set new gain_itr
		gain_itr = gain_bucket[partition].insert(pair <int, CELL*> (gain, this));		
	}
};


//reads the cell file "design.nodes" and fills the array "cells" with all the 
//instances in the design. Also finds their size based on the type of instance.
void readCellsFile()
{
	printf("\nUsing design in %s\n", file_path);

	FILE *cells_file;
	char cells_file_path[100];
	strcpy(cells_file_path, file_path);
	cells_file = fopen(strcat(cells_file_path, "design.nodes"), "r");	
	if (cells_file == NULL) { fprintf(stderr, "Can't open cells_file!\n");	exit(1); }
		
	printf("Reading cells file...\n");
			
	//read through cell instance file design.nodes
	while(!feof(cells_file))
	{
		char inst_name [100];
		char cell_type[20];
		fscanf(cells_file, "%s %s", &inst_name, &cell_type);	 
		if(feof(cells_file)) break;								
		
		//cell_areas map used to find areas with the cell name as keys
		map <string, int> cell_areas;
		cell_areas["FDRE"] = 5;
		cell_areas["LUT6"] = 7;
		cell_areas["LUT5"] = 6;
		cell_areas["LUT4"] = 5;
		cell_areas["LUT3"] = 4;
		cell_areas["LUT2"] = 3;
		cell_areas["LUT1"] = 2;
		cell_areas["CARRY8"] = 34;
		cell_areas["DSP48E2"] = 429;
		cell_areas["RAMB36E2"] = 379;
		cell_areas["BUFGCE"] = 3;
		cell_areas["IBUF"] = 2;
		cell_areas["OBUF"] = 2; 	
		
		//create a new cell and initialize default values
		CELL* current_cell;
		current_cell = new CELL;
		
		vector <NET*> current_net_list;
					
		strcpy(current_cell->name, inst_name);	
		strcpy(current_cell->type, cell_type);				
		current_cell->area = cell_areas[cell_type];			
		current_cell->gain = 0;
		current_cell->fixed = false;
		current_cell->locked = false;
		current_cell->net_list = current_net_list;		
		current_cell->partition = 0; //place all cells in partition 0 initially
		current_cell->best_partition = 0;
		
		partition_area[0] += current_cell->area;
		
		cells.push_back(current_cell);
		
		cells_map[current_cell->name] = current_cell;	
	}
	
	fclose(cells_file);
	printf("done!\n");	
	return;
}


//reads the nets file "design.nets". Fills in "nets" array with all the nets in the design.
void readNetsFile()
{	
	FILE *nets_file;			
	char cells_file_path[100];
	strcpy(cells_file_path, file_path);
	nets_file = fopen(strcat(cells_file_path, "design.nets"), "r");			
	if (nets_file == NULL) { fprintf(stderr, "Can't open nets_file!\n");	exit(1); }
	
	printf("Reading nets file...\n");	
			
	while(!feof(nets_file))
	{
		//read in the net name and number of pins
		NET* current_net;
		current_net = new NET;
			
		char net_name[100];
		fscanf(nets_file, "%*s %s %d", &net_name, &current_net->num_pins);
		if(feof(nets_file)) break;
		
		strcpy(current_net->name, net_name);
		current_net->cutstate = false; //assume not cut
		current_net->critical = false;									
	
		//read each cell in the net, add to net_list and cell_list per FM paper
		for(int c = 0; c < current_net->num_pins; c++)
		{
			char cell_name [100];
			char pin_name [100];	

			fscanf(nets_file, "%s %s", &cell_name, &pin_name);
						
			//avoid adding duplicates
			if(cells_map[cell_name]->net_list.empty() ||
			   cells_map[cell_name]->net_list.back() != current_net )
			{
				//add net to cell's net_list.				
				cells_map[cell_name]->net_list.push_back(current_net);

				//add the new cell to the net's cell_list
				current_net->cell_list.push_back(cells_map[cell_name]);	
			}					
		}
		
		nets.push_back(current_net);			
		fscanf(nets_file, "%*s"); //scan in the "endnet" line				
	}

	fclose(nets_file);
	printf("done!\n");
	return;			
}


void initialPartition()
{
	srand(time(0)); //use current time for random seed				
	int index;
	
	//randomly move cells into other partition until approximately half area in each
	while((float)partition_area[1] / (float)(partition_area[0] + partition_area[1]) < .5 )
	{
		//select a random cell to move to partition 1
		index = rand() % cells.size();
		if(cells[index]->partition == 0)
		{	
			cells[index]->partition = 1;
			cells[index]->best_partition = 1;
			partition_area[0] -= cells[index]->area;
			partition_area[1] += cells[index]->area;		
		}
	}
	
	printf("\nInitial Partition:\n");		
	printPartitionAreas();		
}

void printPartitionAreas()
{
	//print partition information
	printf("Partition 0 area: %d\n", partition_area[0]);
	printf("Partition 1 area: %d\n", partition_area[1]);
	printf("Area Ratio: %.3f/%.3f = %.3f\n\n", (float)partition_area[0] / (float)(partition_area[0] + partition_area[1]), (float)partition_area[1] / (float)(partition_area[0] + partition_area[1]), (float)partition_area[0] / (float)(partition_area[1]));
}


//computes the cutsize of the partitions
//stores in the global variable "cut_size"
void findCutsize()
{
	cut_size = 0;	
	
	//loop through each net. If there is at least one cell in each partition, add 1 to cutsize	
	for(int i = 0; i < nets.size(); i++)
	{		
		nets[i]->partition_count[0] = 0;
		nets[i]->partition_count[1] = 0;
		
		for(int j = 0; j < nets[i]->cell_list.size(); j++)
		{
			//count how many cells in each partition on this net
			nets[i]->partition_count[nets[i]->cell_list[j]->partition]++;			
		}
		
		if(nets[i]->partition_count[0] > 0 && nets[i]->partition_count[1] > 0)
			cut_size++;
	}				
	printf("Cutsize: %d\n", cut_size);
}

int printCutsize()  //finds and returns the cutsize from scratch
{
	int my_cut_size = 0;
	int my_partition_count[2];	
	
	//loop through each net. If there is at least one cell in each partition, add 1 to cutsize	
	for(int i = 0; i < nets.size(); i++)
	{		
		my_partition_count[0] = 0;
		my_partition_count[1] = 0;
		
		for(int j = 0; j < nets[i]->cell_list.size(); j++)
		{
			//count how many cells in each partition on this net
			my_partition_count[nets[i]->cell_list[j]->partition]++;			
		}

		if(my_partition_count[0] && my_partition_count[1])
			my_cut_size++;
	}		
		
	return my_cut_size;
}


//runs through all nets in the design. Computes the cell's potential swap gain based on neighbors
//gain is the decrease in cutsize that would be obtained by swapping a cell
void computeGains()
{
	gain_bucket[0].clear();
	gain_bucket[1].clear();
	
	//look through all cells, compute gain for each cell 
	for(int i = 0; i < cells.size(); i++)
	{	
		cells[i]->gain = 0;
		
			
		//for each net the cell is in, check if there is any gain for that net
		for (int j = 0; j < cells[i]->net_list.size(); j++)
		{	
			NET* current_net = cells[i]->net_list[j];				
						
			//if there is exactly one cell in the "From" partition, increase gain
			if(current_net->partition_count[cells[i]->partition] == 1)
			{	
				cells[i]->gain++;
				current_net->critical = true;
			}
			
			//if there is exactly zero cells in the "To" partition, decrease gain
			if(current_net->partition_count[!cells[i]->partition] == 0)
			{	
				cells[i]->gain--;
				current_net->critical = true;
			}			
		}
		
		//add this cell to a gain bucket, and initialize gain_itr				
		cells[i]->gain_itr = gain_bucket[cells[i]->partition].insert(pair <int, CELL*> (cells[i]->gain,  cells[i]));	
	} 			
}

void saveBestSolution()
{	
	//keep track of the best solution found
	min_cut_size = cut_size;
			
	for(int i = 0; i < cells.size(); i++)
		cells[i]->best_partition = cells[i]->partition;
				
	//printf("SAVED new best solution: (cutsize = %d)\n", cut_size);	
}



//run after each pass
//resets the state of the cells to the best solution found so far
//also resets other attributes to set up for the next pass
void recallBestSolution()
{
	cut_size = min_cut_size;
	
	//reset the cells partitions to the best solution
	//unlock all cells
	for(int i = 0; i < cells.size(); i++)
	{
		cells[i]->partition = cells[i]->best_partition;
		cells[i]->locked = false;		
	}				

	//reset the partition counts for the best solution
	for(int i = 0; i < nets.size(); i++)
	{	
		//reset nets	
		nets[i]->partition_count[0] = 0;
		nets[i]->partition_count[1] = 0;
		
		nets[i]->has_locked[0] = 0; 
		nets[i]->has_locked[1] = 0;
		
		for(int j = 0; j < nets[i]->cell_list.size(); j++)
		{
			//count how many cells in each partition on this net
			nets[i]->partition_count[nets[i]->cell_list[j]->partition]++;			
		}		
	}
	
	//recalculate the gains for the best solution
	computeGains(); 
	
	//set the partition areas to the best solution
	partition_area[0] = 0;
	partition_area[1] = 0;
	
	for(int i = 0; i < cells.size(); i++)
	{		
		partition_area[cells[i]->partition] += cells[i]->area;
	}	
}


//performs 1 pass of the Fidducia-Mathyeses bi-partition algorithm
//assumes all net, cell, gain data has been properly initialized
//returns true if the pass resulted in a better cutsize
//returns false if no improvement was found
bool FMPartitionPass()
{
	printf("\n****BEGIN PARTITION PASS #%d****\n", pass_count);
	
	int starting_cut_size = cut_size;
	int swap_count = 0;
	positive_gains_exhausted = false; //true after the initial gain bucket has been exhausted of positive gains
					//used to prevent updating the best_partition too often	
	int saves_performed = 0;
	
	//perform moves until no legal moves remain
	while(performNextMove())
	{		
		swap_count++;
		
		//if a solution with better cut_size was found AND all positive gains have been exhausted
		if(positive_gains_exhausted && cut_size <= min_cut_size-20)
		{
			saveBestSolution();
			saves_performed++;
		}				
	} 	
		    	
	//revert back to the best solution that was found this pass
	recallBestSolution();
		
	printf("****PASS %d COMPLETE: ", pass_count);
	findCutsize();
	printTimestamp();
	
	if(starting_cut_size == cut_size) //if no improvement was made i.e. stuck in local minimum
		return false;
	else return true;
}


//select the next base_cell
//if moving the cell would imbalance the partitioning, pick another cell
//finally, move the base_cell to the other partition
bool performNextMove()
{	
	if(gain_bucket[0].size() == 0 && gain_bucket[1].size() == 0) //if no more cell in gain bucket, end the pass
		return false;
		
	CELL* base_cell;	
	multimap <int, CELL*>::iterator itr [2]; //iterator for the gain buckets	
	
	if(gain_bucket[0].size() != 0)
		itr[0] = --gain_bucket[0].end(); //pointer to highest gain cell in partition 0
	if(gain_bucket[1].size() != 0)
		itr[1] = --gain_bucket[1].end();
		
	int from_partition;
			
	//if there is nothing left in one of the gain_buckets, use the other bucket
	if(gain_bucket[0].size() == 0) //if no cells in bucket 0
		from_partition = 1;
	else if(gain_bucket[1].size() == 0) //if no cells in bucket 1
		from_partition = 0;
	else from_partition = base_partition;
	
	base_partition = !base_partition; //pick from alternating partitions to make large nets dead faster		
	
	base_cell = itr[from_partition]->second; //pick the highest gain cell as default 

//if the highest gain cell would cause imbalance, pick another cell that doesn't imbalance
if(!cellCanMove(base_cell))
{
	printf("Cannot move %s (Partition: %d, Area: %d) due to imbalance.\tP0 Area: %d\tP1 Area: %d\n", base_cell->name, base_cell->partition, base_cell->area, partition_area[0], partition_area[1]);
		//pick the highest gain cell from the other partition
		from_partition = !from_partition;
		base_cell = itr[from_partition]->second;
		
		if(!cellCanMove(base_cell)) //if the highest gain cell in BOTH partitions would cause an imbalance, no legal moves. End the pass
		{
			printf("Cannot move %s (Partition: %d, Area: %d) due to imbalance.\tP0 Area: %d\tP1 Area: %d\n", base_cell->name, base_cell->partition, base_cell->area, partition_area[0], partition_area[1]);
			printf("No legal moves remain. Both would cause imbalance!\n");
			printf("gain_bucket[0].size() = %d\t gain_bucket[1].size() = %d\n", gain_bucket[0].size(), gain_bucket[1].size());
			return false;
		}
}

	(*base_cell).changePartition();
	
	//return true indicating a successful move
	return true;
}

//check if the indicated cell can legally be moved to the other partition without causing imbalance
bool cellCanMove(CELL* base_cell)
{
//printf("\nCan cell %s move?", base_cell->name);
	if(base_cell == NULL) { printf("base_cell is NULL!\n"); return false;}
	return (float)(partition_area[base_cell->partition] - base_cell->area) / (float)(partition_area[0] + partition_area[1]) >= MIN_CUT_RATIO;
}

//prints how long the program has executed and time since previous timestamp
void printTimestamp()
{

	printf("TIME ELAPSED: %.02f\t(%.02f since previous)\n", (std::clock() - start_time) / (double) CLOCKS_PER_SEC, (std::clock() - timestamp) / (double) CLOCKS_PER_SEC);
	
	timestamp = clock();
}


//printToFile prints output information about the FM partition after completion, 
//including which cells go in which partition for the best cutsize solution found.
void printToFile()
{
	FILE *out_file;
	char out_file_path[100];
	strcpy(out_file_path, file_path);	
	out_file = fopen(strcat(out_file_path, "fmpartition.nodes"), "w");		
			
	printf("\nPrinting results to file: %s\n", out_file_path);
	fprintf(out_file, "FM partition performed on: %s\n\n", file_path);
	fprintf(out_file, "Execution time: %.02f\t\t(%d passes performed)\n", (std::clock() - start_time) / (double) CLOCKS_PER_SEC, (std::clock() - timestamp) / (double) CLOCKS_PER_SEC, pass_count);
	fprintf(out_file, "Starting cut: %d\n", initial_partition_cut_size);
	fprintf(out_file, "Final cut: %d\n", cut_size);
	fprintf(out_file, "Percent change: %.02f\n", (float)(initial_partition_cut_size - cut_size) / (float)initial_partition_cut_size);
	
	//compute ratio cut = cut_size / (cell_count[0] * cell_count[1])
	long cell_count[2] = {0, 0};
	for(int i = 0; i < cells.size(); i++)
		cell_count[cells[i]->partition]++;
		
	fprintf(out_file, "Ratio cut: %.08f\t(%d / (%d*%d))\n", ((float)(cut_size)/(float)(cell_count[0]))/ (float)(cell_count[1]), cut_size, cell_count[0], cell_count[1]);
	
	fprintf(out_file, "\nPartitioned cell list:\n");
	fprintf(out_file, "\tname\t\tcell type (area)\tpartition\n");
	for(int i = 0; i < cells.size(); i++)
	{
		fprintf(out_file, "%15s %15s(%d)\t\t%d\n", cells[i]->name, cells[i]->type, cells[i]->area, cells[i]->best_partition);
	}
}

int main(int argc, char *argv[])
{	
	start_time = std::clock();
	timestamp = std::clock();
	
	//add the example number from the argument to the filepath	
	if(argc > 1) strcat(file_path, argv[1]);
	else strcat(file_path, "1"); //defaults to example 1
	strcat(file_path, "/");
	
	readCellsFile();
	readNetsFile();
	initialPartition();
	
	printf("Initial ");
	findCutsize();
	min_cut_size = cut_size; //initialize min_cut_size
	initial_partition_cut_size = cut_size; //remember the initial cut_size
	
	computeGains();	
	printTimestamp();
	printf("************************************************\n");
		
	int max_num_passes = 20;
	
	//perform passes until local minimum is found, or max_num_passes is exceeded
	while(FMPartitionPass() && pass_count < max_num_passes)
	{
		pass_count++;		
	}		
		
	printf("\n************************************************\nFinal partition areas:\n");
	printPartitionAreas();
	
	printf("Final ");
	findCutsize();
	
	printToFile();	
	
}
