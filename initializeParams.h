#include <vector>

#define PI 3.14159265
int initializeParameters();

// structure to hold coordinates of point
struct point
{
	double x;
	double y;
	double z;
};
// structure to hold referred row and column of a point
struct position
{
	int row;
	int column;
	int depth;
};
// structure to hold data of a growth candidate
struct candidate
{
	int id;
	char branch_type;
	/* valid branch types
	T - terminal
	L - lateral
    t - terminal negative
    l - lateral negative
	*/
};
// structure to hold data related to branch
struct branch
{
	int nodeId;
	point currNode;
	point prevNode;
};

struct branchNew
{
	int nodeId;
	point currNode;
	point parentNode;
	int strahler_number;
	char branch_type;
	/* valid branch types
	T - terminal
	L - lateral
    t - terminal negative
    l - lateral negative
	X - reserved for trunk
	*/

};
// structure to hold data related to nodes in topological tree
struct node
{
	int strahler_number;
	char branch_type;
	/* valid branch types
	T - terminal
	L - lateral
    t - terminal negative
    l - lateral negative
	*/
};

// to be user in future - binary tree data structure
struct node1
{
	node *terminalNode;
	node *lateralNode;
	point coordinates;
	int strahler_number;
	char branch_type;
	/* valid branch types
	T - terminal
	L - lateral
    t - terminal negative
    l - lateral negative
	*/
};

struct obstacle
{
	double obsXStart;
	double obsYStart;
	double obsZStart;

	double obsXEnd;
	double obsYEnd;
	double obsZEnd;
};


namespace constants
{
	// Definitions of the Grid
	extern int VERT_LENGTH; // rows
	extern int HORZ_LENGTH; // cols
	extern int DEPTH_LENGTH; // depth

	// SOR parameters
	extern double CUTOFF;
    extern double ORP; // over reaction parameter

	// tree growth bundry
	extern double MIN_X;
	extern double MAX_X;
	extern double MIN_Y;
	extern double MAX_Y;
	extern double MIN_Z;
	extern double MAX_Z;

	extern double NETA;

	extern int MAX_LEVELS; // total number of branching levels

	//Brancing angles
	extern double THETA;
	extern double BETA;
	extern double DELTA;

	extern double BUD_LENGTH;
	extern double BR_SCALE_FACTOR;	// branch length scaler factor
	extern double BUD_RADIUS;

// Random factors added to growth parameters
	extern double R_THETA;
	extern double R_BETA;
	extern double R_DELTA;
	extern double R_SF; // scale factor

	extern bool SHOW_TREE_GROWTH;

	extern double SOURCE_X_START;
	extern double SOURCE_X_END;
	extern double SOURCE_Y_START;
	extern double SOURCE_Y_END;
	extern double SOURCE_Z_START;
	extern double SOURCE_Z_END;

	extern std::vector<obstacle> OBSTACLES;
};