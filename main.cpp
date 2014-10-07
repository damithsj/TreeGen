#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>
#include <fstream>

// openGL & GLUT includes
#include <GL/glut.h>
#include <GL/gl.h>

#include "definitions.h"
#include "initializeParams.h"
#include "evalSOR.h"
#include "mathUtil.h"
#include "calcSubBranchCoords.h"
#include "openglUtil.h"
#include "generateBranch.h"
#include "gluCylinders.h"

using namespace std;
using namespace constants;

#pragma region Global variables

dVec3D grid;
dVec3D source_distrib;

vector<candidate> growthCandidates; // growth candidates
vector<point> treePoints; // binary tree which holds all points related to tree
vector<branch> treeData; // stores only the grown branches, typically used to draw tree
vector<branchNew> treeDataNew; // stores only the grown branches, typically used to draw tree
vector<node> treeTopology; // stores the structural data related to tree

clock_t esplased_time;

vector<double> brLengthRatio; // branch length ratios (to be used with strahler number)

bool continueDevelopment = true; // break condition for generating tree

int branchCount = 0; // used as part of the file name in writePOVoutputIncremental
#pragma endregion

#pragma region private method declarations
void initializeParams();
void addLightSource();
void addObstacles();
position locatePosition(point p);
int calcParentNodeID(int id);
void generateTree(); // tree generation logic
candidate evalEnvInput();
void addGrowthCandidates(int id);
void addTreePoints(point newNode, int id, char branchType);
void drawTree();
double calcBranchRadius(branch br);
void drawSource();
void drawFloor();
void drawBoundry();
void writePOVoutput();
void writePOVoutputIncremental(); // to be used when the show_tree_growth = true
void drawLeaf(ofstream& outFile, branch br);
void drawTwig(ofstream& outFile, branch br);
void drawTwigLeaf(ofstream& outFile, branch br); // to be used with drawTwig()

bool validateStrahlerNumber(int id);
int updateSubTreeStrahlerNumbers(int id);
int updateParentStrahlerNumber(int id);
bool regenerateSubTree(int startNode); // re generate the sub tree after Strahler Number change
bool isTerminalBranch(branch br); // need for POV output
bool identifyCollision(point point_);
point locateSourceMiddle(double branchLength_);
void drawSourcePOV(); //draw light source for POV

void generateTreeNew();
int getBranchById1(int id);
branchNew getBranchById(int id);
bool validateStrahlerNumberNew(int id);
int updateSubTreeStrahlerNumbersNew(int id);
void drawTreeNew();
double calcBranchRadiusNew(branchNew br);
int updateParentStrahlerNumberNew(int id);
void drawTwigNew(ofstream& outFile, branchNew br);
void drawTwigLeafNew(ofstream& outFile, branchNew br);
void writePOVoutputIncrementalNew();
candidate evalEnvInputNew();
int calcParentNodeIDNew(int id);
#pragma endregion

int main(int argc, char** argv)
{
	int retval;

	retval  =  initializeParameters();
	
	glutInit(&argc, argv);              // initialize glut
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(800,800);
	glutCreateWindow("my Awesome 3D Tree");       // create a window
	initRendering(800,800); //Initialize rendering

	if (SHOW_TREE_GROWTH)
	{
		//glutDisplayFunc(generateTree);
		//generateTree();
		generateTreeNew();
	}
	else
	{
		//generateTree();
		generateTreeNew();
		glutDisplayFunc(drawTree);           // say what the 'display' func is
	}
	glutReshapeFunc(handleResize);
	glutKeyboardFunc(handleKeypress);
	glutSpecialFunc( handleSpecialKeyPress );		// Handles "special" keyboard keys
	//glutFullScreen(); //full screen
	//setColors();

	glutMainLoop();
	return EXIT_SUCCESS;


	/*
	glutInit(&argc, argv);              // initialize glut
	glutInitWindowSize(700,700);	
	glutCreateWindow("my Awesome 2D Tree");       // create a window
	
	if (SHOW_TREE_GROWTH)
	{
		glutDisplayFunc(generateTree); 
	}
	else
	{
		glutDisplayFunc(drawTree);           // say what the 'display' func is
	}
	
	glutReshapeFunc(handleResize);
	glutKeyboardFunc(handleKeypress);

	//setColors();

	glutMainLoop(); */
	return EXIT_SUCCESS;
}


void generateTreeNew()
{
#pragma region private variables	
	branchNew branch_; // branch instance
	position pos_; // position of the branch node
	grid = dVec3D(HORZ_LENGTH, dVec2D( VERT_LENGTH, dVec1D (DEPTH_LENGTH) ));
	source_distrib = dVec3D(HORZ_LENGTH, dVec2D( VERT_LENGTH, dVec1D (DEPTH_LENGTH) )); //resource distribution matrix
	candidate candidate_; // growth candidate

	int parent_;
	branchNew parentBranch_;
	node node_;
	point newNode_;
	bool collisionExist_ = false;
	int subBranchRootId_;

#pragma endregion
	if (continueDevelopment)
	{
		esplased_time = clock();	
		initializeParams(); // initialize global parameters befoe using!!! otherwise they will grow with re drawing in openGL	

		addLightSource();
		addObstacles();
		// evaluate light distribution
		cout << "Evaluating potential distribution..." << endl;
		evalSOR (grid, source_distrib);

		cout << "Start generating tree..." << endl;
		
		// store the trunk points
		branch_.parentNode.x = 0.0f;
		branch_.parentNode.y = 0.0f;
		branch_.parentNode.z = 0.0f;

		branch_.currNode = locateSourceMiddle(BUD_LENGTH);
		branch_.nodeId = 1;
		branch_.strahler_number = 1; // assign 1 every time when adding a node. 
		branch_.branch_type = 'X'; // trunk does not have a branch type

		treeDataNew.push_back(branch_);	

		//update the grid with node point and re calculate the potential
		pos_ = locatePosition(branch_.currNode);
		grid[pos_.column][pos_.row][pos_.depth] = 0;
		source_distrib[pos_.column][pos_.row][pos_.depth] = 1;

		evalSOR (grid, source_distrib);

		// add growth candidates
		candidate_.id = 2;
		candidate_.branch_type = 'T';	
		growthCandidates.push_back(candidate_);

		candidate_.id = 3;
		candidate_.branch_type = 'l';
		growthCandidates.push_back(candidate_);

		// add the trunk in the output
		writePOVoutputIncrementalNew();
		branchCount += 1;

		while (continueDevelopment)
		{
			candidate_ = evalEnvInputNew(); // evaluate enviroment conditions and choose a candidate to grow
			if (candidate_.id >= pow(2.0, MAX_LEVELS)) continueDevelopment = false; // safety break condition

			// find the parent branch id of the candidate		
			parent_ = calcParentNodeIDNew(candidate_.id);
			// get the corresponding data for parent branch
			parentBranch_ = getBranchById(parent_);

			// strahler number should be 1 for every new branch
			node_.branch_type = candidate_.branch_type;
			node_.strahler_number = 1;

			// generate coordinates of new branch
			newNode_ = generateBranch(parentBranch_.parentNode, parentBranch_.currNode, node_);			

			

			// check whether new branch collides with obstacle
			// if there's a collision, we cannot grow current candidate into a branch
			collisionExist_ = identifyCollision(newNode_);
			if (!collisionExist_)
			{
				// add branch to tree
				branch_.nodeId = candidate_.id;
				branch_.parentNode = parentBranch_.currNode;
				branch_.currNode = newNode_;
				branch_.strahler_number = 1;
				branch_.branch_type = candidate_.branch_type;
				treeDataNew.push_back(branch_);
				//update the grid with node point and re calculate the potential
				pos_ = locatePosition(newNode_);
				grid[pos_.column][pos_.row][pos_.depth] = 0;
				source_distrib[pos_.column][pos_.row][pos_.depth] = 1;
				evalSOR (grid, source_distrib);

				// add growth candidates
				addGrowthCandidates(candidate_.id);

				if (validateStrahlerNumberNew(candidate_.id) ==  true) // too bad. inner branch properties needs to change 
				{
					// update strahler_number of parent nodes until root is reached			
				subBranchRootId_ = updateSubTreeStrahlerNumbersNew(candidate_.id);
				// re generate tree starting from above ID
				// check collision exist
				}
				if (SHOW_TREE_GROWTH)
				{
					drawTreeNew();
					writePOVoutputIncrementalNew();
					branchCount += 1;
				}

			}
			else
			{
				continue;
			}
		}
		}
		else
	{
		drawTreeNew();
	}

	//writePOVoutput();
}

void generateTree()
{
#pragma region private variables
	// coordinates of newly generated branch
	//float subX;
	//float subY;

	int iterator = 0;// iterator	
	point point_; // point in the grid
	position pos_; // position of the branch node
	candidate candidate_; // growth candidate

	int parent_;
	point startPointParent_;
	point endPointParent_;
	point newNode_;

	branch branch_;
	node node_;
	int subBranchRootId_;

	grid = dVec3D(HORZ_LENGTH, dVec2D( VERT_LENGTH, dVec1D (DEPTH_LENGTH) ));
	source_distrib = dVec3D(HORZ_LENGTH, dVec2D( VERT_LENGTH, dVec1D (DEPTH_LENGTH) )); //resource distribution matrix

	//positions of the trunk
	float xA = 0.0f;
	float yA = 0.0f;
	float zA = 0.0f;


	point pointB_;

	brLengthRatio = linspace(1.0, MAX_LEVELS, MAX_LEVELS);
	int currLevel_ = 2;

	bool collisionExist_ = false;

	// backups
	dVec3D grid_bkp;
	dVec3D source_distrib_bkp;
	vector<candidate> growthCandidates_bkp;
	vector<point> treePoints_bkp;
	vector<branch> treeData_bkp;
	vector<node> treeTopology_bkp;
#pragma endregion
	if (continueDevelopment)
	{
		esplased_time = clock();	
		initializeParams(); // initialize global parameters befoe using!!! otherwise they will grow with re drawing in openGL	

		addLightSource();
		addObstacles();
		// evaluate light distribution
		cout << "Evaluating potential distribution..." << endl;
		evalSOR (grid, source_distrib);

		cout << "Start generating tree..." << endl;

		// store the branch points
		point_.x = xA;
		point_.y = yA;
		point_.z = zA;
		treePoints.push_back(point_);

		// add dummy node first to tree topology
		node_.strahler_number = 0;
		node_.branch_type = 'X';
		treeTopology.push_back(node_);

		// add the first branch to tree, continue after assigning point_ to B
		branch_.nodeId = 1;
		branch_.prevNode = point_;

		pointB_ = locateSourceMiddle(BUD_LENGTH);

		point_.x = pointB_.x;
		point_.y = pointB_.y;
		point_.z = pointB_.z;
		treePoints.push_back(point_);

		// add the topological data
		node_.strahler_number = 1; // assign 1 every time when adding a node. 
		node_.branch_type = 'X';
		// this will be changed later when node becomes a inner node
		treeTopology.push_back(node_);

		// adding the branch to tree continues...
		branch_.currNode = point_;	
		treeData.push_back(branch_);

		//update the grid with node point and re calculate the potential
		pos_ = locatePosition(point_);

		grid[pos_.column][pos_.row][pos_.depth] = 0;
		source_distrib[pos_.column][pos_.row][pos_.depth] = 1;

		evalSOR (grid, source_distrib);

		// add growth candidates
		candidate_.id = 2;
		candidate_.branch_type = 'T';	
		growthCandidates.push_back(candidate_);

		candidate_.id = 3;
		candidate_.branch_type = 'l';
		growthCandidates.push_back(candidate_);

		// add the trunk in the output
		writePOVoutputIncremental();
		branchCount += 1;

		while (continueDevelopment)
		{		

			//continueDevelopment = false; // set this temperarily
			candidate_ = evalEnvInput(); // evaluate enviroment conditions and choose a candidate to grow
			//cout << "chosen candidate id: " << candidate_.id << endl;

			// ok now we know which branch to grow. why waiting...

			//if (candidate_.id >= pow(2.0, MAX_LEVELS-1)) continueDevelopment = false; // safety break condition
			if (candidate_.id >= pow(2.0, MAX_LEVELS)) continueDevelopment = false; // safety break condition
			// calculate coordinates of the new node

			// find the parent branch id of the candidate		
			parent_ = calcParentNodeID(candidate_.id);		

			// find the start and end points of parent branch

			startPointParent_ = treePoints.at(calcParentNodeID(parent_));

			endPointParent_ = treePoints.at(parent_);

			// strahler number should be 1 for every new branch
			node_.branch_type = candidate_.branch_type;
			node_.strahler_number = 1;

			// generate coordinates of new branch
			newNode_ = generateBranch(startPointParent_, endPointParent_, node_);

			// add branch to tree
			branch_.nodeId = candidate_.id;
			branch_.prevNode = endPointParent_;
			branch_.currNode = newNode_;

			// need to check whether newly added branch will go through the obstacle
			// or if the tree re grneration will go through the obstacle.
			// in tht case we have to giveup this candidate

			//take backups of public variables if we want to restore them if new branch 
			//collide with obstacle
			grid_bkp = grid;
			source_distrib_bkp = source_distrib;

			growthCandidates_bkp =  growthCandidates; // growth candidates
			treePoints_bkp = treePoints; // binary tree which holds all points related to tree
			treeData_bkp = treeData; // stores only the grown branches, typically used to draw tree
			treeTopology_bkp = treeTopology; // stores the structural data related to tree	

			treeData.push_back(branch_);
			// add new node to binary tree
			//treePoints.push_back(newNode_);
			addTreePoints(newNode_, candidate_.id, candidate_.branch_type);
			// TODO
			//Okay, now we added the branch. it's time to calculate new stahler number
			//things need to done...
			//look at the sibiling
			//if strahler_number of sibiling > strahler_number currr node
			//do nothing
			//else (==)
			// loop till find root
			//traverse to parent
			// assign new strahler_number
			// then calculate positions of every node according to their strahler_number

			if (validateStrahlerNumber(candidate_.id) ==  false) // ok falks, we don't need to worry about strahler_number
			{
				// check collision exist
				collisionExist_ = identifyCollision(newNode_);
				if (!collisionExist_)
				{
					//update the grid with node point and re calculate the potential
					pos_ = locatePosition(newNode_);
					grid[pos_.column][pos_.row][pos_.depth] = 0;
					source_distrib[pos_.column][pos_.row][pos_.depth] = 1;
					evalSOR (grid, source_distrib);
				}
			}
			else // too bad
			{
				//cout << "growth candidate ID: " << candidate_.id << endl;
				// update strahler_number of parent nodes until root is reached			
				subBranchRootId_ = updateSubTreeStrahlerNumbers(candidate_.id);

				//cout << "tree starting from " << subBranchRootId_ << " needs to re calculate" << endl;

				// re generate tree starting from above ID
				// check collision exist
				collisionExist_ = regenerateSubTree(subBranchRootId_);
				if (!collisionExist_)
				{
					evalSOR (grid, source_distrib);
				}
			}
			// if collision exist, do not continue with the development of the selected bud.
			// restore the public variables with backups and continue
			if (collisionExist_)
			{
				grid = grid_bkp;
				source_distrib = source_distrib_bkp;
				//growthCandidates = growthCandidates_bkp;  no pint of restoring growth candidates since it will add the
														  //current bud again and cause the algorithm to consider it again
				treePoints = treePoints_bkp;
				treeData = treeData_bkp;
				treeTopology = treeTopology_bkp;
				continue;
			}		

			if (SHOW_TREE_GROWTH)
			{
				drawTree();
				writePOVoutputIncremental();
				branchCount += 1;
			}

			// write the output in order to prevent crash (2013-08-03)
			//cout << treeTopology.size() << endl;
			//if (treeTopology.size() == pow(2.0, currLevel_))
			//{
			//	currLevel_ = currLevel_ + 1;
			//
			//	if (currLevel_ >= MAX_LEVELS - 5)
			//	{
			//		writePOVoutput();
			//	}
			//}
		}

		esplased_time = clock() - esplased_time;
		//cout << "Total time taken to generate tree: " << (float)(esplased_time/CLOCKS_PER_SEC) << " seconds" << endl;
		printf ("Total time taken to generate tree: %f seconds.\n",(float)(esplased_time)/CLOCKS_PER_SEC);
	}
	else
	{
		drawTree();
	}

	writePOVoutput();

}


int getBranchById1(int id)
{
	for (int i = 0; i < treeDataNew.size(); i++)
	{
		if (treeDataNew.at(i).nodeId == id)
		{
			return i;
			break;
		}
	}
	//return treeDataNew[0];
}

branchNew getBranchById(int id)
{
	point dummyPoint_;
	dummyPoint_.x = 0.0;
	dummyPoint_.y = 0.0;
	dummyPoint_.z = 0.0;

	branchNew dummyBranch_;
	dummyBranch_.nodeId = id;
	dummyBranch_.currNode = dummyPoint_;
	dummyBranch_.parentNode = dummyPoint_;
	dummyBranch_.strahler_number = 0;
	dummyBranch_.branch_type = 'X';


	for (int i = 0; i < treeDataNew.size(); i++)
	{
		//cout << "i= " << i << endl;
		//cout << "treedata size= " << treeDataNew.size() << endl;
		//cout << "treeDataNew.at(i).nodeId= " << treeDataNew.at(i).nodeId << endl;
		//cout << "id= " << id << endl;
		//cout << "-------------------------------" << endl;
		if (treeDataNew.at(i).nodeId == id)
		{
			return treeDataNew.at(i);
			break;
		}
	}
	return dummyBranch_;
}

int getRowIdbyNodeId(int id)
{
	for (int i = 0; i < treeDataNew.size(); i++)
	{
		if (treeDataNew[i].nodeId == id)
		{
			return i;
		}
	}
}

int getTreeSize()
{
	int maxNodeId = 0;
	for (int i = 0; i < treeDataNew.size(); i++)
	{
		if (treeDataNew[i].nodeId > maxNodeId)
		{
			maxNodeId = treeDataNew[i].nodeId;
		}
	}
	return maxNodeId;
}
// initialize parameters to zero
void initializeParams()
{
	// initialize arrays to zero;
	for (int i=0; i < HORZ_LENGTH; i++) 
	{
		for (int j=0; j < VERT_LENGTH; j++) 
		{
			for (int k=0; k < DEPTH_LENGTH; k++) 
			{
				grid[i][j][k] = 0.0;
				source_distrib[i][j][k] = 0.0;
			}
		}
	}
	growthCandidates.clear();
	treePoints.clear();
	treeDataNew.clear();
	seed(); // initialize random number generator
	drawSourcePOV(); // generate the POV file containig the source
}

void addLightSource()
{
	//convert the light source configurations to positions in grid
	point point__;
	position startPosition__;
	position endPosition__;

	point__.x = SOURCE_X_START;
	point__.y = SOURCE_Y_START;
	point__.z = SOURCE_Z_START;

	startPosition__ = locatePosition(point__);

	point__.x = SOURCE_X_END;
	point__.y = SOURCE_Y_END;
	point__.z = SOURCE_Z_END;

	endPosition__ = locatePosition(point__);

	//cout << "start col: " << startPosition__.column << endl;
	//cout << "start row: " << startPosition__.row << endl;
	//cout << "start depth: " << startPosition__.depth << endl;

	//cout << "end col: " << endPosition__.column << endl;
	//cout << "end row: " << endPosition__.row << endl;
	//cout << "end depth: " << endPosition__.depth << endl;

	// put the light source
	for (int i = startPosition__.column; i <= endPosition__.column; i++)
	{
		for (int j = startPosition__.row; j <= endPosition__.row; j++)
		{
			for (int k = startPosition__.depth; k >= endPosition__.depth; k--)
			{
				grid[i][j][k] = 1;
				source_distrib[i][j][k] = 1;
			}			
		}		
	}	
}

void addObstacles()
{
//convert the light source configurations to positions in grid

	for(int cnt = 0; cnt < OBSTACLES.size(); cnt++)
	{
	point point__;
	position startPosition__;
	position endPosition__;
	obstacle obs_;

	obs_ = OBSTACLES.at(cnt);
	point__.x = obs_.obsXStart;
	point__.y = obs_.obsYStart;
	point__.z = obs_.obsZStart;
	startPosition__ = locatePosition(point__);

	point__.x = obs_.obsXEnd;
	point__.y = obs_.obsYEnd;
	point__.z = obs_.obsZEnd;
	endPosition__ = locatePosition(point__);

	for (int i = startPosition__.column; i <= endPosition__.column; i++)
	{
		for (int j = startPosition__.row; j <= endPosition__.row; j++)
		{
			for (int k = startPosition__.depth; k >= endPosition__.depth; k--)
			{
				grid[i][j][k] = 0;
				source_distrib[i][j][k] = 0;
			}			
		}		
	}	
	}



	/*point point__;
	position startPosition__;
	position endPosition__;

	point__.x = OBS_X_START;
	point__.y = OBS_Y_START;
	point__.z = OBS_Z_START;

	startPosition__ = locatePosition(point__);

	point__.x = OBS_X_END;
	point__.y = OBS_Y_END;
	point__.z = OBS_Z_END;

	endPosition__ = locatePosition(point__);*/

	//cout << "start col: " << startPosition__.column << endl;
	//cout << "start row: " << startPosition__.row << endl;
	//cout << "start depth: " << startPosition__.depth << endl;

	//cout << "end col: " << endPosition__.column << endl;
	//cout << "end row: " << endPosition__.row << endl;
	//cout << "end depth: " << endPosition__.depth << endl;

	// put the obstacle
	/*for (int i = startPosition__.column; i <= endPosition__.column; i++)
	{
		for (int j = startPosition__.row; j <= endPosition__.row; j++)
		{
			for (int k = startPosition__.depth; k >= endPosition__.depth; k--)
			{
				grid[i][j][k] = 0;
				source_distrib[i][j][k] = 0;
			}			
		}		
	}	*/
}

//return position (row, col) of a given point
position locatePosition(point p)
{
	vector<double> mapX = linspace(MIN_X,MAX_X,HORZ_LENGTH);
	vector<double> mapY = linspace(MIN_Y,MAX_Y,VERT_LENGTH);
	vector<double> mapZ = linspace(MAX_Z,MIN_Z,DEPTH_LENGTH);
	position pos;

	int pointRow = 0;
	int pointCol = 0;
	int pointDepth = DEPTH_LENGTH - 1;

	while (mapX.at(pointCol) < p.x)
	{
		pointCol = pointCol + 1;
	}
	while (mapY.at(pointRow) < p.y)
	{
		pointRow = pointRow + 1;
	}
	while (mapZ.at(pointDepth)  < p.z)
	{
		pointDepth = pointDepth - 1;
	}
	pos.column = pointCol;
	pos.row = pointRow;
	pos.depth = pointDepth;

	return pos;
}

int calcParentNodeID(int id)
{
	if ( id % 2 == 0) // even candidate
		return id/2;
	else
		return (id - 1)/2;
}

int calcParentNodeIDNew(int id)
{
	if ( id % 2 == 0) // even candidate
		return (id)/2;
	else
		return (id - 1)/2;
}

// this function evaluates the growth candidates 
// against the enviromental input and decide
// which candidate to grow into a branch
candidate evalEnvInput()
{
#pragma region private variables
	int parent_;
	point startPoint_;
	point endPoint_;
	point predictNode_;
	position predictPos_;
	candidate returnVal;

	struct candidatePot_ { candidate cVal; double potential; }; // structure to hold candidate data with potential
	vector<candidatePot_> growthCandidatesPot_; //vector to hold all candidate data with potentials

	double potAcc_ = 0; // accumalated potential
	double potTemp_; // temp variable to hold potential value
	candidatePot_ candPot__;
	double probAcc_ = 0;// accumalated probability of occurance

	branchNew branchN_;
#pragma endregion

	for (int i = 0; i < growthCandidates.size(); i++) // repeat for each candidate
	{
		// find the parent branch id of the candidate
		parent_ = calcParentNodeID(growthCandidates.at(i).id);
		//cout << "parent_: " << parent_  << endl;
		// find the parent branch and allocate start and end points of parent branch
		branchN_ = treeDataNew.at(0);
		startPoint_ = branchN_.parentNode;
		endPoint_ = branchN_.currNode;

		//if (parent_ % 2 == 0) // even candidate
		//	startPoint_ = treePoints.at(parent_/2);
		//else
		//	startPoint_ = treePoints.at((parent_ - 1)/2);

		//endPoint_ = treePoints.at(parent_);

		// predict  the projection of development candidate
		predictNode_ = generateBranch(startPoint_, endPoint_, growthCandidates.at(i).branch_type, brLengthRatio.at(0)); // use the length ratio as 1.0. 
		// we assume all buds equal to lenght of terminal branch																										

		// find the predicted node position in grid
		predictPos_ =  locatePosition(predictNode_);

		// get potential from grid
		potTemp_ = grid[predictPos_.column][predictPos_.row][predictPos_.depth];

		//cout << "row: " << predictPos_.row << " col: " << predictPos_.column << endl;

		// store the potential value
		//candPot__.cVal = growthCandidates.at(i); // candidate  *(growthCandidates.begin() + i)
		candPot__.cVal = *(growthCandidates.begin() + i);
		candPot__.potential = potTemp_; // potential
		growthCandidatesPot_.push_back(candPot__);
		//cout << "pot temp: "<<potTemp_ << endl;
		// accumalate the poential
		potAcc_ += pow(abs(potTemp_), NETA); // use absolute value to avoid the NaN error

		//cout << "pot acc: "<<potAcc_ << endl;
	}

	for (int i = 0; i < growthCandidatesPot_.size(); i++)
	{
		if (growthCandidatesPot_.at(i).potential != 0) // do not grow the branc if the corresponding potential is zero (eg: inside obstacle)
		{
			probAcc_ += pow(abs(growthCandidatesPot_.at(i).potential), NETA)/potAcc_;
			//cout << growthCandidatesPot_.at(i).cVal.id << endl;
			//cout << growthCandidatesPot_.at(i).potential << endl;
			//cout << potAcc_ << endl;
			// choose candidate to grow using Monte-Carlo
			if (evalProbability(probAcc_)) // condition to grow satisfied
			{
				returnVal = growthCandidates.at(i);			
				growthCandidates.erase (growthCandidates.begin() + i); // remove chosen node from candidates

				//return growthCandidatesPot_.at(i).cVal; // return the candidate
				return returnVal;
				break;
			}
		}

	}
}

candidate evalEnvInputNew()
{
#pragma region private variables
	int parent_;
	point startPoint_;
	point endPoint_;
	point predictNode_;
	position predictPos_;
	candidate returnVal;
	struct candidatePot_ { candidate cVal; double potential; }; // structure to hold candidate data with potential
	vector<candidatePot_> growthCandidatesPot_; //vector to hold all candidate data with potentials
	double potAcc_ = 0; // accumalated potential
	double potTemp_; // temp variable to hold potential value
	candidatePot_ candPot__;
	double probAcc_ = 0;// accumalated probability of occurance
	branchNew branch_;
	branchNew parentBranch_;
#pragma endregion
	
	for (int i = 0; i < growthCandidates.size(); i++) // repeat for each candidate
	{
		// find the parent branch id of the candidate
		if ( growthCandidates.at(i).id % 2 == 0) // even candidate
			parent_ = (growthCandidates.at(i).id)/2;
		else
			parent_ = (growthCandidates.at(i).id - 1)/2;

		// find the start and end points of parent branch
		parentBranch_  = getBranchById(parent_);

		startPoint_ = parentBranch_.parentNode;
		endPoint_ = parentBranch_.currNode;

		
		//cout << "startPoint_ " << startPoint_.x << " " << startPoint_.y << " " << startPoint_.z << endl;
		//cout << "endPoint_ " << endPoint_.x << " " << endPoint_.y << " " << endPoint_.z << endl;

		//cout << "growthCandidates.at(i) " << growthCandidates.at(i).id << " " << growthCandidates.at(i).branch_type << endl;
		// predict  the projection of development candidate
		predictNode_ = generateBranch(startPoint_, endPoint_, growthCandidates.at(i).branch_type, 1.0); // use the length ratio as 1.0. 
		// we assume all buds equal to lenght of terminal branch
		
		//cout << "predictNode_ " << predictNode_.x << " " << predictNode_.y << " " << predictNode_.z << endl;
		//cout << "-------------------------------" << endl;
		// find the predicted node position in grid
		predictPos_ =  locatePosition(predictNode_);

		// get potential from grid
		potTemp_ = grid[predictPos_.column][predictPos_.row][predictPos_.depth];

		//cout << "row: " << predictPos_.row << " col: " << predictPos_.column << endl;

		// store the potential value
		//candPot__.cVal = growthCandidates.at(i); // candidate  *(growthCandidates.begin() + i)
		candPot__.cVal = *(growthCandidates.begin() + i);
		candPot__.potential = potTemp_; // potential
		growthCandidatesPot_.push_back(candPot__);
		//cout << "pot temp: "<<potTemp_ << endl;
		// accumalate the poential
		potAcc_ += pow(abs(potTemp_), NETA); // use absolute value to avoid the NaN error

		//cout << "pot acc: "<<potAcc_ << endl;
	}

	for (int i = 0; i < growthCandidatesPot_.size(); i++)
	{
		if (growthCandidatesPot_.at(i).potential != 0) // do not grow the branc if the corresponding potential is zero (eg: inside obstacle)
		{
			probAcc_ += pow(abs(growthCandidatesPot_.at(i).potential), NETA)/potAcc_;
			//cout << growthCandidatesPot_.at(i).cVal.id << endl;
			//cout << growthCandidatesPot_.at(i).potential << endl;
			//cout << potAcc_ << endl;
			// choose candidate to grow using Monte-Carlo
			if (evalProbability(probAcc_)) // condition to grow satisfied
			{
				returnVal = growthCandidates.at(i);			
				growthCandidates.erase (growthCandidates.begin() + i); // remove chosen node from candidates

				//return growthCandidatesPot_.at(i).cVal; // return the candidate
				return returnVal;
				break;
			}
		}

	}
}

void addGrowthCandidates(int id)
{
	candidate candidate_; // growth candidate

	if (id % 2 == 0) // even node
	{
		candidate_.id = 2*id;
		candidate_.branch_type = 'L';	
		growthCandidates.push_back(candidate_);

		candidate_.id = 2*id + 1;
		candidate_.branch_type = 't';
		growthCandidates.push_back(candidate_);
	}
	else
	{
		candidate_.id = 2*id;
		candidate_.branch_type = 'T';	
		growthCandidates.push_back(candidate_);

		candidate_.id = 2*id + 1;
		candidate_.branch_type = 'l';
		growthCandidates.push_back(candidate_);
	}
}

void addTreePoints(point newNode, int id, char branchType)
{
	point dummyPoint_;
	dummyPoint_.x = 0.0;
	dummyPoint_.y = 0.0;
	dummyPoint_.z = 0.0;

	node dummyNode_;
	dummyNode_.strahler_number = 0; // assign strahler number 0 for dormant nodes
	dummyNode_.branch_type = 'X'; // assign strahler number 0 for dormant nodes

	if (id >= treePoints.size())	treePoints.resize(id + 1, dummyPoint_); // resize the tree if growth candidate is greater

	if (id >= treeTopology.size())	treeTopology.resize(id + 1, dummyNode_); // resize the tree topology if growth candidate is greater/
	// treeTopology should be a mirror of treePoints

	dummyNode_.strahler_number = 1; // assign strahler number 1 for new branch
	dummyNode_.branch_type = branchType; // assign strahler number 1 for new branch
	treePoints[id] = newNode;
	treeTopology[id] = dummyNode_;
}

void drawTree()
{
	double radius = BUD_RADIUS;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear last draw
	glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
	glLoadIdentity();

	glTranslatef(movementX, movementY, movementZ);
	glRotatef( rotateAngle, 0.0, 1.0, 0.0 ); 
	glRotatef( azimuthAngle, 1.0, 0.0, 0.0 );

	for (int i = 0;  i < treeData.size(); i++)
	{
		// calculate radius for the branch
		radius = calcBranchRadius(treeData.at(i));
		//generate branch
		renderCylinder_convenient(treeData.at(i).prevNode.x, treeData.at(i).prevNode.y, treeData.at(i).prevNode.z, treeData.at(i).currNode.x, treeData.at(i).currNode.y, treeData.at(i).currNode.z,radius,16);

		//cout << "start x: " << treeData.at(i).prevNode.x << endl;
		//cout << "start y: " << treeData.at(i).prevNode.y << endl;
		//cout << "start z: " << treeData.at(i).prevNode.z << endl;

		//cout << "end x: " << treeData.at(i).currNode.x << endl;
		//cout << "end y: " << treeData.at(i).currNode.y << endl;
		//cout << "end z: " << treeData.at(i).currNode.z << endl;

		//cout << "radius: " << radius << endl;

		//cout << "----------------------------------------------------------" << endl;		
	}
	
	drawSource();
	drawFloor();
	drawBoundry();

	glutSwapBuffers(); //Send the 3D scene to the screen

	//cout << "Total time taken to draw tree: " << esplased_time/CLOCKS_PER_SEC << " seconds" << endl;
	//printf ("Total time taken to draw tree: %f seconds.\n",(float)(esplased_time)/CLOCKS_PER_SEC);
}

void drawTreeNew()
{
	double radius = BUD_RADIUS;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // clear last draw
	glMatrixMode(GL_MODELVIEW); //Switch to the drawing perspective
	glLoadIdentity();

	glTranslatef(movementX, movementY, movementZ);
	glRotatef( rotateAngle, 0.0, 1.0, 0.0 ); 
	glRotatef( azimuthAngle, 1.0, 0.0, 0.0 );

	for (int i = 0;  i < treeDataNew.size(); i++)
	{
		// calculate radius for the branch
		radius = calcBranchRadiusNew(treeDataNew.at(i));
		//generate branch
		renderCylinder_convenient(treeDataNew.at(i).parentNode.x, treeDataNew.at(i).parentNode.y,treeDataNew.at(i).parentNode.z, treeDataNew.at(i).currNode.x, treeDataNew.at(i).currNode.y, treeDataNew.at(i).currNode.z,radius,16);

		//cout << "start x: " << treeData.at(i).prevNode.x << endl;
		//cout << "start y: " << treeData.at(i).prevNode.y << endl;
		//cout << "start z: " << treeData.at(i).prevNode.z << endl;

		//cout << "end x: " << treeData.at(i).currNode.x << endl;
		//cout << "end y: " << treeData.at(i).currNode.y << endl;
		//cout << "end z: " << treeData.at(i).currNode.z << endl;

		//cout << "radius: " << radius << endl;

		//cout << "----------------------------------------------------------" << endl;		
	}
	
	drawSource();
	drawFloor();
	drawBoundry();

	glutSwapBuffers(); //Send the 3D scene to the screen

	//cout << "Total time taken to draw tree: " << esplased_time/CLOCKS_PER_SEC << " seconds" << endl;
	//printf ("Total time taken to draw tree: %f seconds.\n",(float)(esplased_time)/CLOCKS_PER_SEC);
}

double calcBranchRadius(branch br)
{
	node node_;
	int branch_sth_no_;

	node_ = treeTopology.at(br.nodeId);
	branch_sth_no_ = node_.strahler_number;

	// changed according to the formula given by Viennot. et.al (ACM 1989)
	//return BUD_RADIUS*branch_sth_no_;
	
	double c2 = 0.67; // constant
	double alpha_ = 1.15;
	double beta_ = 1.5;

	//return BUD_RADIUS* c2 * pow(branch_sth_no_, alpha_); //  for polynomial growth (added bud radius)
	//return c2 * pow(beta_, branch_sth_no_); // for exponential growth
	return BUD_RADIUS;
}
double calcBranchRadiusNew(branchNew br)
{
	node node_;
	int branch_sth_no_;

	//node_ = treeTopology.at(br.nodeId);
	branch_sth_no_ = br.strahler_number;

	// changed according to the formula given by Viennot. et.al (ACM 1989)
	//return BUD_RADIUS*branch_sth_no_;
	
	double c2 = 0.67; // constant
	double alpha_ = 1.15;
	//double beta_ = 1.5;
	double beta_ = 1.1;

	return BUD_RADIUS* c2 * pow(branch_sth_no_, alpha_); //  for polynomial growth (added bud radius)
	//return c2 * pow(beta_, branch_sth_no_); // for exponential growth
	//return BUD_RADIUS;
}

// this function returns true if sibiling
// Strahler Number is equal to current Strahler Number
bool validateStrahlerNumber(int id)
{

	// whatever node generated, look at the sibiling and parent to decide update others or not

	int sibiling_strahler_nr;
	int current_strahler_nr;
	int parent_strahler_nr;

	current_strahler_nr = treeTopology.at(id).strahler_number;
	parent_strahler_nr = treeTopology.at(calcParentNodeID(id)).strahler_number;


	/*current node is a even node and it's the final node. therefore no sibilings. should not update the sub tree*/
	if ((id % 2 == 0) && (id + 1 == treeTopology.size()))  
	{
		return false;
	}

	// calculate sibiling strahler number
	if (id % 2 == 0) sibiling_strahler_nr = treeTopology.at(id + 1).strahler_number;
	else sibiling_strahler_nr = treeTopology.at(id - 1).strahler_number;

	//cout << "---sth validation: current id " << id << endl;
	//cout << "---sth validation: current_strahler_nr " << current_strahler_nr << endl;
	//cout << "---sth validation: parent_strahler_nr " << parent_strahler_nr << endl;
	//cout << "---sth validation: sibiling_strahler_nr " << sibiling_strahler_nr << endl;

	if ((current_strahler_nr == sibiling_strahler_nr) && (current_strahler_nr == parent_strahler_nr))
	{
		return true;
	}

	if ((current_strahler_nr > sibiling_strahler_nr) && (current_strahler_nr > parent_strahler_nr))
	{
		return true;
	}

	if ((current_strahler_nr == sibiling_strahler_nr) && (current_strahler_nr > parent_strahler_nr))
	{
		return true;
	}
	// terminal node
	//if ((current_strahler_nr == sibiling_strahler_nr) && (current_strahler_nr ==  1) && (current_strahler_nr >= parent_strahler_nr)) 
	//{
	//	return true;
	//}

	//if ()

	//// intermidiate node
	//if ((current_strahler_nr >= sibiling_strahler_nr) && (current_strahler_nr >= parent_strahler_nr)) 
	//{		
	//	return true;
	//}

	return false;
}

bool validateStrahlerNumberNew(int id)
{

	// whatever node generated, look at the sibiling and parent to decide update others or not

	int sibiling_strahler_nr;
	int current_strahler_nr;
	int parent_strahler_nr;
	branchNew branch_;
	branchNew parentBranch_;

	branch_ = getBranchById(id);
	parentBranch_ = getBranchById(calcParentNodeIDNew(id));

	current_strahler_nr = branch_.strahler_number;
	parent_strahler_nr = parentBranch_.strahler_number;


	/*current node is a even node and it's the final node. therefore no sibilings. should not update the sub tree*/
	if ((id % 2 == 0) && (id + 1 == getTreeSize()))  
	{
		return false;
	}

	// calculate sibiling strahler number
	if (id % 2 == 0) sibiling_strahler_nr = getBranchById(id + 1).strahler_number;
	else sibiling_strahler_nr = getBranchById(id - 1).strahler_number;

	//cout << "---sth validation: current id " << id << endl;
	//cout << "---sth validation: current_strahler_nr " << current_strahler_nr << endl;
	//cout << "---sth validation: parent_strahler_nr " << parent_strahler_nr << endl;
	//cout << "---sth validation: sibiling_strahler_nr " << sibiling_strahler_nr << endl;

	if ((current_strahler_nr == sibiling_strahler_nr) && (current_strahler_nr == parent_strahler_nr))
	{
		return true;
	}

	if ((current_strahler_nr > sibiling_strahler_nr) && (current_strahler_nr > parent_strahler_nr))
	{
		return true;
	}

	if ((current_strahler_nr == sibiling_strahler_nr) && (current_strahler_nr > parent_strahler_nr))
	{
		return true;
	}
	// terminal node
	//if ((current_strahler_nr == sibiling_strahler_nr) && (current_strahler_nr ==  1) && (current_strahler_nr >= parent_strahler_nr)) 
	//{
	//	return true;
	//}

	//if ()

	//// intermidiate node
	//if ((current_strahler_nr >= sibiling_strahler_nr) && (current_strahler_nr >= parent_strahler_nr)) 
	//{		
	//	return true;
	//}

	return false;
}

int updateSubTreeStrahlerNumbers(int id)
{	
	while ((id > 1) && (validateStrahlerNumber(id)))
	{
		//cout << "=====update subtree of "<< id << endl;
		id = updateParentStrahlerNumber(id);		
	}

	return id;
}

int updateSubTreeStrahlerNumbersNew(int id)
{	
	while ((id > 1) && (validateStrahlerNumberNew(id)))
	{
		//cout << "=====update subtree of "<< id << endl;
		id = updateParentStrahlerNumberNew(id);		
	}

	return id;
}

// returns updated parent ID
int updateParentStrahlerNumber(int id)
{
	int nodeToUpdate_ = calcParentNodeID(id); // here we update the parent's sth_no
	int sibiling_strahler_nr;
	int current_strahler_nr;
	int nodeToUpdate_strahler_nr;
	node node_ = treeTopology.at(nodeToUpdate_);
	
	// get current node strahler number
	current_strahler_nr = treeTopology.at(id).strahler_number;

	// get the strahler number of node to update;
	nodeToUpdate_strahler_nr = treeTopology.at(nodeToUpdate_).strahler_number;

	// get sibiling strahler number
	if (id % 2 == 0) sibiling_strahler_nr = treeTopology.at(id + 1).strahler_number;
	else sibiling_strahler_nr = treeTopology.at(id - 1).strahler_number;


	//cout << "+++++++++++++++current_strahler_nr: " << current_strahler_nr << endl;
	//cout << "+++++++++++++++sibiling_strahler_nr: " << sibiling_strahler_nr << endl;
	//cout << "+++++++++++++++nodeToUpdate_strahler_nr: " << nodeToUpdate_strahler_nr << endl;

	if (current_strahler_nr == sibiling_strahler_nr)
	{
		//cout << "--------------previous strahler_number of " << nodeToUpdate_ << " th node "<< node_.strahler_number << endl;
		node_.strahler_number = treeTopology.at(id).strahler_number + 1; //increase the strahler number by one
		treeTopology[nodeToUpdate_] = node_;
		
		//cout << "--------------new strahler_number of " << nodeToUpdate_ << " th node "<< node_.strahler_number << endl;
	}
	else if (current_strahler_nr > nodeToUpdate_strahler_nr)		
	{
		//cout << "+++++++++++++++previous strahler_number: " << node_.strahler_number << endl;
		node_.strahler_number = current_strahler_nr; //assign the strahler number of current node
		treeTopology[nodeToUpdate_] = node_;
		//cout << "+++++++++++++++new strahler_number: " << node_.strahler_number << endl;
	}

	
	//cout << "new strahler_number: " << node_.strahler_number << endl;
	return nodeToUpdate_;
	
}

// returns updated parent ID
int updateParentStrahlerNumberNew(int id)
{
	int nodeToUpdate_ = calcParentNodeID(id); // here we update the parent's sth_no
	int sibiling_strahler_nr;
	int current_strahler_nr;
	int nodeToUpdate_strahler_nr;
	//node node_ = treeTopology.at(nodeToUpdate_);
	branchNew branch_= getBranchById(nodeToUpdate_);
	
	// get current node strahler number
	current_strahler_nr = getBranchById(id).strahler_number;

	// get the strahler number of node to update;
	nodeToUpdate_strahler_nr = getBranchById(nodeToUpdate_).strahler_number;

	// get sibiling strahler number
	if (id % 2 == 0) sibiling_strahler_nr = getBranchById(id + 1).strahler_number;
	else sibiling_strahler_nr = getBranchById(id - 1).strahler_number;


	//cout << "+++++++++++++++current_strahler_nr: " << current_strahler_nr << endl;
	//cout << "+++++++++++++++sibiling_strahler_nr: " << sibiling_strahler_nr << endl;
	//cout << "+++++++++++++++nodeToUpdate_strahler_nr: " << nodeToUpdate_strahler_nr << endl;

	if (current_strahler_nr == sibiling_strahler_nr)
	{
		//cout << "--------------previous strahler_number of " << nodeToUpdate_ << " th node "<< node_.strahler_number << endl;
		//node_.strahler_number = treeTopology.at(id).strahler_number + 1; //increase the strahler number by one
		branch_.strahler_number = getBranchById(id).strahler_number + 1; //increase the strahler number by one
		
		treeDataNew[getRowIdbyNodeId(nodeToUpdate_)] = branch_;
		
		//cout << "--------------new strahler_number of " << nodeToUpdate_ << " th node "<< node_.strahler_number << endl;
	}
	else if (current_strahler_nr > nodeToUpdate_strahler_nr)		
	{
		//cout << "+++++++++++++++previous strahler_number: " << node_.strahler_number << endl;
		branch_.strahler_number = current_strahler_nr; //assign the strahler number of current node
		treeDataNew[getRowIdbyNodeId(nodeToUpdate_)] = branch_;
		//cout << "+++++++++++++++new strahler_number: " << node_.strahler_number << endl;
	}

	
	//cout << "new strahler_number: " << node_.strahler_number << endl;
	return nodeToUpdate_;
	
}
bool regenerateSubTree(int startNode)
{
	int parent_;
	point dummyPoint_;
	point startPointParent_;
	point endPointParent_;
	point currPointBkp_;
	point updatedPoint_;
	branch branch_;
	position pos_;
	dummyPoint_.x = 0.0;
	dummyPoint_.y = 0.0;
	dummyPoint_.z = 0.0;
	double c1 = 1.6;
	bool collisionExist_ = false;

	for (int i = startNode; i < treePoints.size();  i++)
	{
		if (!((treePoints.at(i).x == dummyPoint_.x) && (treePoints.at(i).y == dummyPoint_.y) && (treePoints.at(i).z == dummyPoint_.z)))
		{

			if (i == 1)
			{
				//cout << "updating trunk should done differently" << endl;
				cout << "root sth no: " << treeTopology.at(i).strahler_number << endl;
				endPointParent_ = treePoints.at(0);
				updatedPoint_ = treePoints.at(i);
				// changed according to the formula given by Viennot. et.al (ACM 1989)
				//updatedPoint_.z = BUD_LENGTH*treeTopology.at(i).strahler_number*BR_SCALE_FACTOR;  
				//updatedPoint_.z = BUD_LENGTH*c1*treeTopology.at(i).strahler_number*BR_SCALE_FACTOR; // for linear growth
				//updatedPoint_ = locateSourceMiddle(BUD_LENGTH*c1*treeTopology.at(i).strahler_number*BR_SCALE_FACTOR); // for linear growth
				//updatedPoint_ = locateSourceMiddle(BUD_LENGTH);
				updatedPoint_ = locateSourceMiddle(BUD_LENGTH*c1*BR_SCALE_FACTOR);
				//updatedPoint_.z = BUD_LENGTH*c1*treeTopology.at(i).strahler_number*treeTopology.at(i).strahler_number*BR_SCALE_FACTOR; // for quadratic growth
					//continue;
			}
			else
			{
				// find the parent branch id of the candidate		
				parent_ = calcParentNodeID(i);
				startPointParent_ = treePoints.at(calcParentNodeID(parent_));
				endPointParent_ = treePoints.at(parent_);

				updatedPoint_ = generateBranchNoRand(startPointParent_, endPointParent_, treeTopology.at(i));
			}
			// check whether updated point reside inside obstacle
			collisionExist_ = identifyCollision(updatedPoint_);
			if (collisionExist_) return collisionExist_;
		
		// take a backup
		currPointBkp_ = treePoints.at(i);

		//revert old value in source distribution
		pos_ = locatePosition(currPointBkp_);
		source_distrib[pos_.column][pos_.row][pos_.depth] = 0;

		// update treePoints
		treePoints[i] = updatedPoint_;
		//cout << "sth no of " << i << " th branch :"<< treeTopology.at(i).strahler_number << endl;
		//cout << "node " << i << " updated" << endl;

		//update the grid with node point and re calculate the potential
		//cout << "updating source distribution should done differently" << endl;
		pos_ = locatePosition(updatedPoint_);
		grid[pos_.column][pos_.row][pos_.depth] = 0;
		source_distrib[pos_.column][pos_.row][pos_.depth] = 1;

		// search and update the tree data also.

		for (int j = 0; j < treeData.size(); j ++)
		{			
			if (treeData.at(j).nodeId == i) // found the node need to change
			{
				branch_ = treeData.at(j);

				branch_.prevNode = endPointParent_;
				branch_.currNode = updatedPoint_;
				treeData[j] = branch_;
				break;
			}
		}	
		}
	}
	return collisionExist_;
}


bool isTerminalBranch(branch br)
{
	node node_;
	int child1_sth_no;
	int child2_sth_no;

	node_ = treeTopology.at(2*br.nodeId); // check child 1
	child1_sth_no = node_.strahler_number;

	node_ = treeTopology.at(2*br.nodeId + 1); // check child 1
	child2_sth_no = node_.strahler_number;

	if (child1_sth_no == 0 && child2_sth_no == 0)
	{
		return true;
	}
	else
	{
		return false;
	}

}

// this method will identify whether given point
// is inside the obstacle
bool identifyCollision(point point_)
{
	bool collition_detected_ = false;
	for(int cnt = 0; cnt < OBSTACLES.size(); cnt++)
	{
		if (((point_.x >= OBSTACLES.at(cnt).obsXStart) && (point_.x <= OBSTACLES.at(cnt).obsXEnd)) && ((point_.y >= OBSTACLES.at(cnt).obsYStart) && (point_.y <= OBSTACLES.at(cnt).obsYEnd)) && ((point_.z >= OBSTACLES.at(cnt).obsZStart) && (point_.z <= OBSTACLES.at(cnt).obsZEnd)))
	{
		return true;
	}
	else
	{
		collition_detected_ = false;
	}

	}

	return collition_detected_;	

}

point locateSourceMiddle(double branchLength_)
{
	double sourceXMiddle_ = SOURCE_X_START + (SOURCE_X_END - SOURCE_X_START)/2;
	double sourceYMiddle_ = SOURCE_Y_START +(SOURCE_Y_END - SOURCE_Y_START)/2;
	double sourceZMiddle_ = SOURCE_Z_START + (SOURCE_Z_END - SOURCE_Z_START)/2;
	point output_;

	double lengthSource_ = sqrt(sourceXMiddle_*sourceXMiddle_ + sourceYMiddle_*sourceYMiddle_ + sourceZMiddle_*sourceZMiddle_);
	
	output_.x = (branchLength_/lengthSource_)*sourceXMiddle_;
	output_.y = (branchLength_/lengthSource_)*sourceYMiddle_;
	output_.z = (branchLength_/lengthSource_)*sourceZMiddle_;
	return output_;
}

void drawFloor()
{
	glBegin(GL_LINES);
	//for(GLfloat x = MIN_X; x <= MAX_X; x+= (MAX_X - MIN_X)/HORZ_LENGTH)
	for(GLfloat x = MIN_X; x <= MAX_X; x+= 0.5) 
	{
		glVertex3f(x, MIN_Y, 0.0f);
		glVertex3f(x, MAX_Y, 0.0f);
	}
	glEnd();

	glBegin(GL_LINES);
	//for(GLfloat y = MIN_Y; y <= MAX_Y; y+= (MAX_Y - MIN_Y)/VERT_LENGTH)
	for(GLfloat y = MIN_Y; y <= MAX_Y; y+= 0.5)
	{
		glVertex3f(MIN_X,  y, 0.0f);
		glVertex3f(MAX_X, y, 0.0f);
	}
	glEnd(); 
}

void drawSource()
{
	glBegin(GL_POLYGON); 
	glVertex3f( SOURCE_X_START, SOURCE_Y_START, SOURCE_Z_START);       // P1
	glVertex3f( SOURCE_X_START, SOURCE_Y_START, SOURCE_Z_END);       // P2
	glVertex3f( SOURCE_X_END, SOURCE_Y_START, SOURCE_Z_END);       // P3
	glVertex3f( SOURCE_X_END, SOURCE_Y_START, SOURCE_Z_START);       // P4 
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f( SOURCE_X_END, SOURCE_Y_START, SOURCE_Z_START);       // P4 
	glVertex3f( SOURCE_X_END, SOURCE_Y_START, SOURCE_Z_END);       // P3
	glVertex3f( SOURCE_X_END, SOURCE_Y_END, SOURCE_Z_END);       // P5
	glVertex3f( SOURCE_X_END, SOURCE_Y_END, SOURCE_Z_START);       // P6
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f( SOURCE_X_END, SOURCE_Y_END, SOURCE_Z_START);       // P6 
	glVertex3f( SOURCE_X_START, SOURCE_Y_END, SOURCE_Z_START);       // P7
	glVertex3f( SOURCE_X_START, SOURCE_Y_END, SOURCE_Z_END);       // P8
	glVertex3f( SOURCE_X_END, SOURCE_Y_END, SOURCE_Z_END);       // P5
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f( SOURCE_X_START, SOURCE_Y_START, SOURCE_Z_START);       // P1
	glVertex3f( SOURCE_X_START, SOURCE_Y_START, SOURCE_Z_END);       // P2
	glVertex3f( SOURCE_X_START, SOURCE_Y_END, SOURCE_Z_END);       // P8
	glVertex3f( SOURCE_X_START, SOURCE_Y_END, SOURCE_Z_START);       // P7
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f( SOURCE_X_START, SOURCE_Y_START, SOURCE_Z_START);       // P1
	glVertex3f( SOURCE_X_START, SOURCE_Y_END, SOURCE_Z_START);       // P7
	glVertex3f( SOURCE_X_END, SOURCE_Y_END, SOURCE_Z_START);       // P6 
	glVertex3f( SOURCE_X_END, SOURCE_Y_START, SOURCE_Z_START);       // P4 
	glEnd();

	glBegin(GL_POLYGON);
	glVertex3f( SOURCE_X_END, SOURCE_Y_START, SOURCE_Z_END);       // P3
	glVertex3f( SOURCE_X_START, SOURCE_Y_START, SOURCE_Z_END);       // P2
	glVertex3f( SOURCE_X_START, SOURCE_Y_END, SOURCE_Z_END);       // P8
	glVertex3f( SOURCE_X_END, SOURCE_Y_END, SOURCE_Z_END);       // P5
	glEnd();
}

void drawBoundry()
{
	// bottom plane
	glBegin(GL_LINES);
	glVertex3f(MIN_X, MIN_Y, MIN_Z);
	glVertex3f(MIN_X, MAX_Y, MIN_Z);

	glVertex3f(MIN_X, MAX_Y, MIN_Z);
	glVertex3f(MAX_X, MAX_Y, MIN_Z);

	glVertex3f(MAX_X, MAX_Y, MIN_Z);
	glVertex3f(MAX_X, MIN_Y, MIN_Z);

	glVertex3f(MAX_X, MIN_Y, MIN_Z);
	glVertex3f(MIN_X, MIN_Y, MIN_Z);
	glEnd();

	// top plane
	glBegin(GL_LINES);
	glVertex3f(MIN_X, MIN_Y, MAX_Z);
	glVertex3f(MIN_X, MAX_Y, MAX_Z);

	glVertex3f(MIN_X, MAX_Y, MAX_Z);
	glVertex3f(MAX_X, MAX_Y, MAX_Z);

	glVertex3f(MAX_X, MAX_Y, MAX_Z);
	glVertex3f(MAX_X, MIN_Y, MAX_Z);

	glVertex3f(MAX_X, MIN_Y, MAX_Z);
	glVertex3f(MIN_X, MIN_Y, MAX_Z);
	glEnd();

	// verticle lines
	glBegin(GL_LINES);
	glVertex3f(MIN_X, MIN_Y, MIN_Z);
	glVertex3f(MIN_X, MIN_Y, MAX_Z);

	glVertex3f(MIN_X, MAX_Y, MIN_Z);
	glVertex3f(MIN_X, MAX_Y, MAX_Z);

	glVertex3f(MAX_X, MAX_Y, MIN_Z);
	glVertex3f(MAX_X, MAX_Y, MAX_Z);

	glVertex3f(MAX_X, MIN_Y, MIN_Z);
	glVertex3f(MAX_X, MIN_Y, MAX_Z);
	glEnd();
}

void writePOVoutput()
{
	double radius;
	string fileName;
	time_t now = time(NULL);
	char timestamp[20];
	strftime(timestamp, 20, "%Y-%m-%d %H-%M-%S", localtime(&now));
	
	fileName.append("../output/").append("tree ").append(timestamp).append(".inc");	
	ofstream outFile;
	outFile.open (fileName);
	//
	//append includes
	outFile << "#include \"colors.inc\"" << endl;	
	outFile << "#include " << "\"include_light_source" << ".inc\"" << endl;
	outFile << "#include " << "\"include_camera" << ".inc\"" << endl;
	outFile << "#include " << "\"include_pov_light_config" << ".inc\"" << endl;
	// Draw tree
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//---------------------------Branch segments-------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;

	outFile << "blob {" << endl;
	outFile << "      threshold 0.07" << endl;


	for (int i = 0;  i < treeData.size(); i++)
	{
		
		//else
		//{
			//draw branches
			radius = calcBranchRadius(treeData.at(i));
			outFile << "      cylinder {" << endl;
			outFile << "      <" << treeData.at(i).prevNode.x << ", " << treeData.at(i).prevNode.y << ", " << treeData.at(i).prevNode.z << ">" << "," << endl;// base point
			outFile << "      <" << treeData.at(i).currNode.x << ", " << treeData.at(i).currNode.y << ", " << treeData.at(i).currNode.z << ">" << "," << endl;// cap point
			outFile << "       " << radius << endl; // radius
			//outFile << "open" << endl; // Remove end caps
			outFile << "       1" << endl; // Blob strength
			//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
			outFile << "       }" << endl;
			outFile << "       sphere {" << endl;
			outFile << "       <" << treeData.at(i).currNode.x << ", " << treeData.at(i).currNode.y << ", " << treeData.at(i).currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
		//}
	}
	outFile << "pigment{color rgb<150/356,81/356,0>} " << endl; //brown color
	outFile << "}" << endl;

	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//--------------------------------Leaves-----------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	for (int i = 0;  i < treeData.size(); i++)
	{
		if (treeData.at(i).nodeId >= pow(2.0, MAX_LEVELS-1))// definite terminal branch 
		{

			//object{branch translate<-1, 0, 0> scale 0.4 translate<0, 0, 2>}
			//outFile << "object { branch " << "translate<-1, 0, 0> scale " << 1*BUD_LENGTH << "translate " << "<" << treeData.at(i).prevNode.x << ", " << treeData.at(i).prevNode.y << ", " << treeData.at(i).prevNode.z << ">" << "}" << endl;
			drawLeaf(outFile, treeData.at(i));
		}
		else if (isTerminalBranch(treeData.at(i)))
		{
			//outFile << "object { branch " << "translate<-1, 0, 0> scale " << 1*BUD_LENGTH << "translate " << "<" << treeData.at(i).prevNode.x << ", " << treeData.at(i).prevNode.y << ", " << treeData.at(i).prevNode.z << ">" << "}" << endl;
			drawLeaf(outFile, treeData.at(i));
		}
	}
}
void writePOVoutputIncremental()
{
	double radius;
	string fileName;
	string fileNamePrefix;
	bool drawTree_ = false;
	char branchCount_[200];
	int sth_no_;

	if (branchCount < 10) 
	{
		fileNamePrefix = "tree_0000000";
		drawTree_ = true; 
	}		
	else if (branchCount < 100) 
	{
		fileNamePrefix = "tree_000000";
		drawTree_ = true;
	}		
	else if (branchCount < 1000) 
	{
		fileNamePrefix = "tree_00000";
		if (branchCount % 5 == 0)
		{
			drawTree_ = true; // draw only 5th tree
		}
	}
	else if (branchCount < 5000) 
	{
		fileNamePrefix = "tree_00000";
		if (branchCount % 10 == 0)
		{
			drawTree_ = true; // draw only 10th tree
		}
	}
	else if (branchCount < 10000) 
	{
		fileNamePrefix = "tree_0000";
		if (branchCount % 100 == 0)
		{
			drawTree_ = true; // draw only 100th tree
		}
	}		
	else if (branchCount < 100000) 
	{
		fileNamePrefix = "tree_000";
		if (branchCount % 1000 == 0)
		{
			drawTree_ = true; // draw only 1000th tree
		}
	}
	else if (branchCount < 1000000) 
	{
		fileNamePrefix = "tree_00";
		if (branchCount % 10000 == 0)
		{
			drawTree_ = true; // draw only 10000th tree
		}
	}
	else if (branchCount < 10000000) 
	{
		fileNamePrefix = "tree_0";
		if (branchCount % 100000 == 0)
		{
			drawTree_ = true; // draw only 100000th tree
		}
	}
	else 
	{
		fileNamePrefix = "tree_";
		if (branchCount % 1000000 == 0)
		{
			drawTree_ = true; // draw only 1000000th tree
		}
	}	
	if (drawTree_)
	{
	itoa(branchCount, branchCount_, 10);
	
	fileName.append("../output/").append(fileNamePrefix).append(branchCount_).append(".inc");	
	ofstream outFile;
	outFile.open (fileName);
	
	//append includes
	outFile << "#include \"colors.inc\"" << endl;

	// Draw tree
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//---------------------------Branch segments-------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;

	outFile << "blob {" << endl;
	outFile << "      threshold 0.07" << endl;


	for (int i = 0;  i < treeData.size(); i++)
	{
		
		//else
		//{
			//draw branches
			radius = calcBranchRadius(treeData.at(i));
			outFile << "      cylinder {" << endl;
			outFile << "      <" << treeData.at(i).prevNode.x << ", " << treeData.at(i).prevNode.y << ", " << treeData.at(i).prevNode.z << ">" << "," << endl;// base point
			outFile << "      <" << treeData.at(i).currNode.x << ", " << treeData.at(i).currNode.y << ", " << treeData.at(i).currNode.z << ">" << "," << endl;// cap point
			outFile << "       " << radius << endl; // radius
			//outFile << "open" << endl; // Remove end caps
			outFile << "       1" << endl; // Blob strength
			//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
			outFile << "       }" << endl;
			outFile << "       sphere {" << endl;
			outFile << "       <" << treeData.at(i).currNode.x << ", " << treeData.at(i).currNode.y << ", " << treeData.at(i).currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
		//}
	}
	outFile << "pigment{color rgb<150/356,81/356,0>} " << endl; //brown color
	outFile << "}" << endl;

	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//--------------------------------Leaves-----------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	for (int i = 0;  i < treeData.size(); i++)
	{
		sth_no_ = treeTopology.at(treeData.at(i).nodeId).strahler_number; // fetch sth no of the branch
		if (sth_no_ == 1)
		{
			drawTwig(outFile, treeData.at(i));
		}		
	}

	outFile.close();
	//cout << fileName.c_str() << endl;

	// create a .pov file containing the .inc file.
	fileName = "";
	fileName.append("../output/").append(fileNamePrefix).append(branchCount_).append(".pov");	
	//ofstream outFile;
	outFile.open (fileName);

	outFile << "#include " << "\"include_light_source" << ".inc\"" << endl;
	outFile << "#include " << "\"include_camera" << ".inc\"" << endl;
	outFile << "#include " << "\"include_pov_light_config" << ".inc\"" << endl;

	outFile << "#include " << "\"" << (char*)fileNamePrefix.c_str() << branchCount_ << ".inc\"" << endl;
	outFile.close();
	}

}
void writePOVoutputIncrementalNew()
{
	double radius;
	string fileName;
	string fileNamePrefix;
	bool drawTree_ = false;
	char branchCount_[200];
	int sth_no_;

	if (branchCount < 10) 
	{
		fileNamePrefix = "tree_0000000";
		drawTree_ = true; 
	}		
	else if (branchCount < 100) 
	{
		fileNamePrefix = "tree_000000";
		drawTree_ = true;
	}		
	else if (branchCount < 1000) 
	{
		fileNamePrefix = "tree_00000";
		if (branchCount % 5 == 0)
		{
			drawTree_ = true; // draw only 5th tree
		}
	}
	else if (branchCount < 5000) 
	{
		fileNamePrefix = "tree_00000";
		if (branchCount % 10 == 0)
		{
			drawTree_ = true; // draw only 10th tree
		}
	}
	else if (branchCount < 10000) 
	{
		fileNamePrefix = "tree_0000";
		if (branchCount % 100 == 0)
		{
			drawTree_ = true; // draw only 100th tree
		}
	}		
	else if (branchCount < 100000) 
	{
		fileNamePrefix = "tree_000";
		if (branchCount % 1000 == 0)
		{
			drawTree_ = true; // draw only 1000th tree
		}
	}
	else if (branchCount < 1000000) 
	{
		fileNamePrefix = "tree_00";
		if (branchCount % 10000 == 0)
		{
			drawTree_ = true; // draw only 10000th tree
		}
	}
	else if (branchCount < 10000000) 
	{
		fileNamePrefix = "tree_0";
		if (branchCount % 100000 == 0)
		{
			drawTree_ = true; // draw only 100000th tree
		}
	}
	else 
	{
		fileNamePrefix = "tree_";
		if (branchCount % 1000000 == 0)
		{
			drawTree_ = true; // draw only 1000000th tree
		}
	}	
	if (drawTree_)
	{
	itoa(branchCount, branchCount_, 10);
	
	fileName.append("../output/").append(fileNamePrefix).append(branchCount_).append(".inc");	
	ofstream outFile;
	outFile.open (fileName);
	
	//append includes
	outFile << "#include \"colors.inc\"" << endl;

	// Draw tree
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//---------------------------Branch segments-------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;

	outFile << "blob {" << endl;
	outFile << "      threshold 0.07" << endl;


	for (int i = 0;  i < treeDataNew.size(); i++)
	{
		
		//else
		//{
			//draw branches
			radius = calcBranchRadiusNew(treeDataNew.at(i));
			outFile << "      cylinder {" << endl;
			outFile << "      <" << treeDataNew.at(i).parentNode.x << ", " << treeDataNew.at(i).parentNode.y << ", " << treeDataNew.at(i).parentNode.z << ">" << "," << endl;// base point
			outFile << "      <" << treeDataNew.at(i).currNode.x << ", " << treeDataNew.at(i).currNode.y << ", " << treeDataNew.at(i).currNode.z << ">" << "," << endl;// cap point
			outFile << "       " << radius << endl; // radius
			//outFile << "open" << endl; // Remove end caps
			outFile << "       1" << endl; // Blob strength
			//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
			outFile << "       }" << endl;
			outFile << "       sphere {" << endl;
			outFile << "       <" << treeDataNew.at(i).currNode.x << ", " << treeDataNew.at(i).currNode.y << ", " << treeDataNew.at(i).currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
		//}
	}
	outFile << "pigment{color rgb<150/356,81/356,0>} " << endl; //brown color
	outFile << "}" << endl;

	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//--------------------------------Leaves-----------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	for (int i = 0;  i < treeDataNew.size(); i++)
	{
		sth_no_ = treeDataNew.at(i).strahler_number; // fetch sth no of the branch
		if (sth_no_ == 1)
		{
			drawTwigNew(outFile, treeDataNew.at(i));
		}		
	}

	outFile.close();
	//cout << fileName.c_str() << endl;

	// create a .pov file containing the .inc file.
	fileName = "";
	fileName.append("../output/").append(fileNamePrefix).append(branchCount_).append(".pov");	
	//ofstream outFile;
	outFile.open (fileName);

	outFile << "#include " << "\"include_light_source" << ".inc\"" << endl;
	outFile << "#include " << "\"include_camera" << ".inc\"" << endl;
	outFile << "#include " << "\"include_pov_light_config" << ".inc\"" << endl;

	outFile << "#include " << "\"" << (char*)fileNamePrefix.c_str() << branchCount_ << ".inc\"" << endl;
	outFile.close();
	}

}

void drawLeaf(ofstream& outFile, branch br)
{
//			(3)
//			/\
//		   /  \
//        /    \
//     (2)\    /(4)
//         \  /
//          \/
//          (1)
//
//
	double xA = br.prevNode.x;
	double yA = br.prevNode.y;
	double zA = br.prevNode.z;
	double xB = br.currNode.x;
	double yB = br.currNode.y;
	double zB = br.currNode.z;
	double subX;
	double subY;
	double subZ;
	
	point point1;
	point point2;
	point point3;
	point point4;

	point1.x = br.currNode.x;
	point1.y = br.currNode.y;
	point1.z = br.currNode.z;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, -30.0*PI/180, 0.0, 0.5);

	point2.x = subX;
	point2.y = subY;
	point2.z = subZ;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 10.0*PI/180, -90.0*PI/180, 1.2);

	point3.x = subX;
	point3.y = subY;
	point3.z = subZ;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 30.0*PI/180, 0.0, 0.5);

	point4.x = subX;
	point4.y = subY;
	point4.z = subZ;

	outFile << "mesh2 {" << endl;
	outFile << "   vertex_vectors {" << endl;
	outFile << "4," << endl;
	outFile << "<" << point1.x << ", " << point1.y << ", " << point1.z << ">" << "," << endl; // 0
	outFile << "<" << point2.x << ", " << point2.y << ", " << point2.z << ">" << "," << endl; // 1
	outFile << "<" << point3.x << ", " << point3.y << ", " << point3.z << ">" << "," << endl; // 2
	outFile << "<" << point4.x << ", " << point4.y << ", " << point4.z << ">" << endl;  // 3
	outFile << "   }" << endl;
	outFile << "   face_indices {" << endl;
	outFile << "2," << endl;
	outFile << "<0, 1, 3>, <1, 2, 3>" << endl;
	outFile << "   }" << endl;
	outFile << "texture { pigment { Green } }" << endl; //color
	outFile << "}" << endl;
}

void drawTwig(ofstream& outFile, branch br)
{
	//	
	//	   @    @<---- leaf
	//		\  /
	//       \/<------ angle: 75, width: bud radius/2
	//       ||
	//       ||<------ branch with sth_no: 1

	branch twig1_;
	branch twig2_;
	double xA = br.prevNode.x;
	double yA = br.prevNode.y;
	double zA = br.prevNode.z;
	double xB = br.currNode.x;
	double yB = br.currNode.y;
	double zB = br.currNode.z;
	double subX;
	double subY;
	double subZ;
	double radius = BUD_RADIUS/4;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, -60.0*PI/180, 0.0, 0.2);
	twig1_.prevNode.x = xB;
	twig1_.prevNode.y = yB;
	twig1_.prevNode.z = zB;

	twig1_.currNode.x = subX;
	twig1_.currNode.y = subY;
	twig1_.currNode.z = subZ;

	// draw twig1
	outFile << "      cylinder {" << endl;
	outFile << "      <" << twig1_.prevNode.x << ", " << twig1_.prevNode.y << ", " << twig1_.prevNode.z << ">" << "," << endl;// base point
	outFile << "      <" << twig1_.currNode.x << ", " << twig1_.currNode.y << ", " << twig1_.currNode.z << ">" << "," << endl;// cap point
	outFile << "       " << radius << endl; // radius
	//outFile << "open" << endl; // Remove end caps
	//outFile << "       1" << endl; // Blob strength
	//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
	outFile << "       }" << endl;
	//outFile << "       sphere {" << endl;
	//outFile << "       <" << twig1_.currNode.x << ", " << twig1_.currNode.y << ", " << twig1_.currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
	
	// draw corresponding leaf
	drawTwigLeaf(outFile, twig1_);

	//------------------------------------------------------------------------------------
	// twig 2
	
	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 60.0*PI/180, 0.0, 0.2);
	twig2_.prevNode.x = xB;
	twig2_.prevNode.y = yB;
	twig2_.prevNode.z = zB;

	twig2_.currNode.x = subX;
	twig2_.currNode.y = subY;
	twig2_.currNode.z = subZ;

	// draw twig2
	outFile << "      cylinder {" << endl;
	outFile << "      <" << twig2_.prevNode.x << ", " << twig2_.prevNode.y << ", " << twig2_.prevNode.z << ">" << "," << endl;// base point
	outFile << "      <" << twig2_.currNode.x << ", " << twig2_.currNode.y << ", " << twig2_.currNode.z << ">" << "," << endl;// cap point
	outFile << "       " << radius << endl; // radius
	//outFile << "open" << endl; // Remove end caps
	//outFile << "       1" << endl; // Blob strength
	//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
	outFile << "       }" << endl;
	//outFile << "       sphere {" << endl;
	//outFile << "       <" << twig2_.currNode.x << ", " << twig2_.currNode.y << ", " << twig2_.currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
	// draw corresponding leaf
	drawTwigLeaf(outFile, twig2_);
}

void drawTwigNew(ofstream& outFile, branchNew br)
{
	//	
	//	   @    @<---- leaf
	//		\  /
	//       \/<------ angle: 75, width: bud radius/2
	//       ||
	//       ||<------ branch with sth_no: 1

	branchNew twig1_;
	branchNew twig2_;
	double xA = br.parentNode.x;
	double yA = br.parentNode.y;
	double zA = br.parentNode.z;
	double xB = br.currNode.x;
	double yB = br.currNode.y;
	double zB = br.currNode.z;
	double subX;
	double subY;
	double subZ;
	double radius = BUD_RADIUS/4;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, -60.0*PI/180, 0.0, 0.2);
	twig1_.parentNode.x = xB;
	twig1_.parentNode.y = yB;
	twig1_.parentNode.z = zB;

	twig1_.currNode.x = subX;
	twig1_.currNode.y = subY;
	twig1_.currNode.z = subZ;

	// draw twig1
	outFile << "      cylinder {" << endl;
	outFile << "      <" << twig1_.parentNode.x << ", " << twig1_.parentNode.y << ", " << twig1_.parentNode.z << ">" << "," << endl;// base point
	outFile << "      <" << twig1_.currNode.x << ", " << twig1_.currNode.y << ", " << twig1_.currNode.z << ">" << "," << endl;// cap point
	outFile << "       " << radius << endl; // radius
	//outFile << "open" << endl; // Remove end caps
	//outFile << "       1" << endl; // Blob strength
	//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
	outFile << "       }" << endl;
	//outFile << "       sphere {" << endl;
	//outFile << "       <" << twig1_.currNode.x << ", " << twig1_.currNode.y << ", " << twig1_.currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
	
	// draw corresponding leaf
	drawTwigLeafNew(outFile, twig1_);

	//------------------------------------------------------------------------------------
	// twig 2
	
	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 60.0*PI/180, 0.0, 0.2);
	twig2_.parentNode.x = xB;
	twig2_.parentNode.y = yB;
	twig2_.parentNode.z = zB;

	twig2_.currNode.x = subX;
	twig2_.currNode.y = subY;
	twig2_.currNode.z = subZ;

	// draw twig2
	outFile << "      cylinder {" << endl;
	outFile << "      <" << twig2_.parentNode.x << ", " << twig2_.parentNode.y << ", " << twig2_.parentNode.z << ">" << "," << endl;// base point
	outFile << "      <" << twig2_.currNode.x << ", " << twig2_.currNode.y << ", " << twig2_.currNode.z << ">" << "," << endl;// cap point
	outFile << "       " << radius << endl; // radius
	//outFile << "open" << endl; // Remove end caps
	//outFile << "       1" << endl; // Blob strength
	//outFile << "texture { pigment{color rgb<150/356,81/356,0>} }" << endl; //brown color
	outFile << "       }" << endl;
	//outFile << "       sphere {" << endl;
	//outFile << "       <" << twig2_.currNode.x << ", " << twig2_.currNode.y << ", " << twig2_.currNode.z << ">" << "," <<  radius << ", -1} " << endl;// //blob sphere with -1 strength
	// draw corresponding leaf
	drawTwigLeafNew(outFile, twig2_);
}

void drawTwigLeaf(ofstream& outFile, branch br)
{
//			(3)
//			/\
//		   /  \
//        /    \
//     (2)\    /(4)
//         \  /
//          \/
//          (1)
//
//
	double xA = br.prevNode.x;
	double yA = br.prevNode.y;
	double zA = br.prevNode.z;
	double xB = br.currNode.x;
	double yB = br.currNode.y;
	double zB = br.currNode.z;
	double subX;
	double subY;
	double subZ;
	
	point point1;
	point point2;
	point point3;
	point point4;

	point1.x = br.currNode.x;
	point1.y = br.currNode.y;
	point1.z = br.currNode.z;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, -30.0*PI/180, 0.0, 0.5/0.2);// 0.2 factor to use the length of branch

	point2.x = subX;
	point2.y = subY;
	point2.z = subZ;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 10.0*PI/180, -90.0*PI/180, 1.2/0.2);

	point3.x = subX;
	point3.y = subY;
	point3.z = subZ;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 30.0*PI/180, 0.0, 0.5/0.2);

	point4.x = subX;
	point4.y = subY;
	point4.z = subZ;

	outFile << "mesh2 {" << endl;
	outFile << "   vertex_vectors {" << endl;
	outFile << "4," << endl;
	outFile << "<" << point1.x << ", " << point1.y << ", " << point1.z << ">" << "," << endl; // 0
	outFile << "<" << point2.x << ", " << point2.y << ", " << point2.z << ">" << "," << endl; // 1
	outFile << "<" << point3.x << ", " << point3.y << ", " << point3.z << ">" << "," << endl; // 2
	outFile << "<" << point4.x << ", " << point4.y << ", " << point4.z << ">" << endl;  // 3
	outFile << "   }" << endl;
	outFile << "   face_indices {" << endl;
	outFile << "2," << endl;
	outFile << "<0, 1, 3>, <1, 2, 3>" << endl;
	outFile << "   }" << endl;
	outFile << "texture { pigment { Green } }" << endl; //color
	outFile << "}" << endl;
}

void drawTwigLeafNew(ofstream& outFile, branchNew br)
{
//			(3)
//			/\
//		   /  \
//        /    \
//     (2)\    /(4)
//         \  /
//          \/
//          (1)
//
//
	double xA = br.parentNode.x;
	double yA = br.parentNode.y;
	double zA = br.parentNode.z;
	double xB = br.currNode.x;
	double yB = br.currNode.y;
	double zB = br.currNode.z;
	double subX;
	double subY;
	double subZ;
	
	point point1;
	point point2;
	point point3;
	point point4;

	point1.x = br.currNode.x;
	point1.y = br.currNode.y;
	point1.z = br.currNode.z;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, -30.0*PI/180, 0.0, 0.5/0.2);// 0.2 factor to use the length of branch

	point2.x = subX;
	point2.y = subY;
	point2.z = subZ;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 10.0*PI/180, -90.0*PI/180, 1.2/0.2);

	point3.x = subX;
	point3.y = subY;
	point3.z = subZ;

	calcSubBranchCoords(xA, yA, zA, xB, yB, zB, subX, subY, subZ, 30.0*PI/180, 0.0, 0.5/0.2);

	point4.x = subX;
	point4.y = subY;
	point4.z = subZ;

	outFile << "mesh2 {" << endl;
	outFile << "   vertex_vectors {" << endl;
	outFile << "4," << endl;
	outFile << "<" << point1.x << ", " << point1.y << ", " << point1.z << ">" << "," << endl; // 0
	outFile << "<" << point2.x << ", " << point2.y << ", " << point2.z << ">" << "," << endl; // 1
	outFile << "<" << point3.x << ", " << point3.y << ", " << point3.z << ">" << "," << endl; // 2
	outFile << "<" << point4.x << ", " << point4.y << ", " << point4.z << ">" << endl;  // 3
	outFile << "   }" << endl;
	outFile << "   face_indices {" << endl;
	outFile << "2," << endl;
	outFile << "<0, 1, 3>, <1, 2, 3>" << endl;
	outFile << "   }" << endl;
	outFile << "texture { pigment { Green } }" << endl; //color
	outFile << "}" << endl;
}

// will be call inside initialize config since we need this file at the begining of rendering
void drawSourcePOV()
{
	string fileName;
	
	fileName.append("../output/").append("include_light_source").append(".inc");	
	ofstream outFile;
	outFile.open (fileName);

	outFile << "#include \"colors.inc\"" << endl;
	// Draw Source
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//---------------------------Light Source----------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "box {" << endl;
	outFile << "<" << SOURCE_X_START << ", " << SOURCE_Y_START << ", " << SOURCE_Z_START << ">" << "," << endl; // Near lower left corner
	outFile << "<" << SOURCE_X_END << ", " << SOURCE_Y_END << ", " << SOURCE_Z_END << ">" << endl; // Far upper right corner
	outFile << "texture { pigment { Yellow } }" << endl; //color
	outFile << "}" << endl;

	// Draw Obstacles
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	outFile << "//------------------------------Obstacles----------------------------------\\" << endl;
	outFile << "//-------------------------------------------------------------------------\\" << endl;
	for(int cnt = 0; cnt < OBSTACLES.size(); cnt++)
	{
		outFile << "box {" << endl;
		outFile << "<" << OBSTACLES.at(cnt).obsXStart << ", " << OBSTACLES.at(cnt).obsYStart << ", " << OBSTACLES.at(cnt).obsZStart << ">" << "," << endl; // Near lower left corner
		outFile << "<" << OBSTACLES.at(cnt).obsXEnd << ", " << OBSTACLES.at(cnt).obsYEnd << ", " << OBSTACLES.at(cnt).obsZEnd << ">" << endl; // Far upper right corner
		outFile << "texture { pigment { Blue } }" << endl; //color
	    outFile << "}" << endl;
	}	

} 
