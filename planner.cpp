
/*=================================================================
 *
 * planner.c
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include <vector>
#include <random>
#include <iostream>
#include <chrono>
#include <unordered_map>

/* Input Arguments */
#define	MAP_IN      prhs[0]
#define	ARMSTART_IN	prhs[1]
#define	ARMGOAL_IN     prhs[2]
#define	PLANNER_ID_IN     prhs[3]

/* Planner Ids */
#define RRT         0
#define RRTCONNECT  1
#define RRTSTAR     2
#define PRM         3

/* Output Arguments */
#define	PLAN_OUT	plhs[0]
#define	PLANLENGTH_OUT	plhs[1]

#define GETMAPINDEX(X, Y, XSIZE, YSIZE) (Y*XSIZE + X)

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#define PI 3.141592654

//the length of each link in the arm (should be the same as the one used in runtest.m)
#define LINKLENGTH_CELLS 10
using namespace std;
typedef struct {
  int X1, Y1;
  int X2, Y2;
  int Increment;
  int UsingYIndex;
  int DeltaX, DeltaY;
  int DTerm;
  int IncrE, IncrNE;
  int XIndex, YIndex;
  int Flipped;
} bresenham_param_t;

const double eps = PI/10.0;
const double step = eps/10.0;
double cost = 0;
const double del = 5.264;
const double gam = 10;

struct node{
	double* angle_;
	node* parent_;
	double g_;
	node(int numofDOFs)
	{
		parent_ = nullptr;
		g_ = 0;
		angle_ = new double[numofDOFs];
	}
	~node()
	{
		parent_ = nullptr;
		delete [] angle_;
	}
};

struct node_star{
	double* angle_;
	double g_;
	node_star* parent_;
	std::vector<node_star*> child_;
	node_star(int numofDOFs) {
		angle_ = new double[numofDOFs];
		parent_ = NULL;
		g_ = 0;
	}
	~node_star()
	{
		parent_ = nullptr;
		delete angle_;
	}
};

struct node_prm{
	double* angle_;
	double g_;
	double h_;
	double f_;
	std::vector<node_prm*> connected;
	node_prm* parent_;
	int comp_index_;
	int node_index_;
	bool visited;
	node_prm(int numofDOFs, int index)
	{
		angle_ = new double[numofDOFs];
		parent_ = nullptr;
		visited = false;
		g_ = 0;
		h_ = 0;
		f_ = g_ + h_;
		comp_index_ = index;
		node_index_ = index;

	}
};

struct compare_func{
	bool operator() (const node_prm* node1, const node_prm* node2) {
		return node1->f_ > node2->f_;
	}
};
void ContXY2Cell(double x, double y, short unsigned int* pX, short unsigned int *pY, int x_size, int y_size)
{
    double cellsize = 1.0;
	//take the nearest cell
	*pX = (int)(x/(double)(cellsize));
	if( x < 0) *pX = 0;
	if( *pX >= x_size) *pX = x_size-1;

	*pY = (int)(y/(double)(cellsize));
	if( y < 0) *pY = 0;
	if( *pY >= y_size) *pY = y_size-1;
}


void get_bresenham_parameters(int p1x, int p1y, int p2x, int p2y, bresenham_param_t *params)
{
  params->UsingYIndex = 0;

  if (fabs((double)(p2y-p1y)/(double)(p2x-p1x)) > 1)
    (params->UsingYIndex)++;

  if (params->UsingYIndex)
    {
      params->Y1=p1x;
      params->X1=p1y;
      params->Y2=p2x;
      params->X2=p2y;
    }
  else
    {
      params->X1=p1x;
      params->Y1=p1y;
      params->X2=p2x;
      params->Y2=p2y;
    }

   if ((p2x - p1x) * (p2y - p1y) < 0)
    {
      params->Flipped = 1;
      params->Y1 = -params->Y1;
      params->Y2 = -params->Y2;
    }
  else
    params->Flipped = 0;

  if (params->X2 > params->X1)
    params->Increment = 1;
  else
    params->Increment = -1;

  params->DeltaX=params->X2-params->X1;
  params->DeltaY=params->Y2-params->Y1;

  params->IncrE=2*params->DeltaY*params->Increment;
  params->IncrNE=2*(params->DeltaY-params->DeltaX)*params->Increment;
  params->DTerm=(2*params->DeltaY-params->DeltaX)*params->Increment;

  params->XIndex = params->X1;
  params->YIndex = params->Y1;
}

void get_current_point(bresenham_param_t *params, int *x, int *y)
{
  if (params->UsingYIndex)
    {
      *y = params->XIndex;
      *x = params->YIndex;
      if (params->Flipped)
        *x = -*x;
    }
  else
    {
      *x = params->XIndex;
      *y = params->YIndex;
      if (params->Flipped)
        *y = -*y;
    }
}

int get_next_point(bresenham_param_t *params)
{
  if (params->XIndex == params->X2)
    {
      return 0;
    }
  params->XIndex += params->Increment;
  if (params->DTerm < 0 || (params->Increment < 0 && params->DTerm <= 0))
    params->DTerm += params->IncrE;
  else
    {
      params->DTerm += params->IncrNE;
      params->YIndex += params->Increment;
    }
  return 1;
}



int IsValidLineSegment(double x0, double y0, double x1, double y1, double*	map,
		   int x_size,
 		   int y_size)

{
	bresenham_param_t params;
	int nX, nY; 
    short unsigned int nX0, nY0, nX1, nY1;

    //printf("checking link <%f %f> to <%f %f>\n", x0,y0,x1,y1);
    
	//make sure the line segment is inside the environment
	if(x0 < 0 || x0 >= x_size ||
		x1 < 0 || x1 >= x_size ||
		y0 < 0 || y0 >= y_size ||
		y1 < 0 || y1 >= y_size)
		return 0;

	ContXY2Cell(x0, y0, &nX0, &nY0, x_size, y_size);
	ContXY2Cell(x1, y1, &nX1, &nY1, x_size, y_size);

    //printf("checking link <%d %d> to <%d %d>\n", nX0,nY0,nX1,nY1);

	//iterate through the points on the segment
	get_bresenham_parameters(nX0, nY0, nX1, nY1, &params);
	do {
		get_current_point(&params, &nX, &nY);
		if(map[GETMAPINDEX(nX,nY,x_size,y_size)] == 1)
            return 0;
	} while (get_next_point(&params));

	return 1;
}

int IsValidArmConfiguration(double* angles, int numofDOFs, double*	map,
		   int x_size, int y_size)
{
    double x0,y0,x1,y1;
    int i;
    
 	//iterate through all the links starting with the base
	x1 = ((double)x_size)/2.0;
    y1 = 0;
	for(i = 0; i < numofDOFs; i++)
	{
		//compute the corresponding line segment
		x0 = x1;
		y0 = y1;
		x1 = x0 + LINKLENGTH_CELLS*cos(2*PI-angles[i]);
		y1 = y0 - LINKLENGTH_CELLS*sin(2*PI-angles[i]);

		//check the validity of the corresponding line segment
		if(!IsValidLineSegment(x0,y0,x1,y1,map,x_size,y_size))
			return 0;
	}
	return 1;    
}
double Distance_between(double* angle1, double* angle2, int numofDOFs) 
{
	double dis;
	double tmp;
	for (int i = 0; i < numofDOFs; i++) {
		tmp = angle1[i] - angle2[i];
		if (tmp > PI) {
			tmp = tmp - 2*PI;
		}
		if (tmp < -PI) {
			tmp = tmp + 2*PI;
		}
		dis = dis + tmp * tmp;
	}
	return sqrt(dis);

}
int isObstacleFree(double* angle1, double* angle2, int numofDOFs, double* map, int x_size, int y_size)
{
	if(!IsValidArmConfiguration(angle2, numofDOFs, map, x_size, y_size)) {
		return 0;
	}
	double delta[5];
	double dis_del = 0;
	double* test_config = new double[numofDOFs];
	for (int i = 0; i < numofDOFs; i++) {
		delta[i] = angle2[i] - angle1[i];
		if (delta[i] < -PI) {
			delta[i] = delta[i] + 2*PI;
		}
		if(delta[i] > PI) {
			delta[i] = delta[i] - 2*PI;
		}
		dis_del = dis_del + delta[i]*delta[i];
	}
	dis_del = sqrt(dis_del);
	double unit_delta[numofDOFs];
	for (int i = 0; i < numofDOFs; i++) {
		unit_delta[i] = delta[i]/dis_del;
	}
	// check obstacle
	for (int i = 0; i < dis_del / step + 1; i++) {
		for (int i = 0; i < numofDOFs; i++) {
			test_config[i] = angle1[i] + i * step * unit_delta[i] ; 
		}
		if(!IsValidArmConfiguration(test_config, numofDOFs, map, x_size, y_size)) {
			return 0;
		}
	}
	return 1;

}

int Extend(std::vector<node*> &tree, node* sample_node, int numofDOFs, double* map, int x_size, int y_size)
{
	// Extend
	// State for Extend
	bool use_near = false;
	bool use_eps = false;
	int extend_status = 0;//  0 for trapped, 1 for reached, 2 for advanced
	// Get the nearest node
	node* nearest_node;
	double dis_s = 0;
	double sum_min = 0;
	for (int i = 0; i < tree.size(); i++) {
		double dis = 0;
		double sum = 0;
		for (int j = 0; j < numofDOFs; j++) {
			dis = tree[i]->angle_[j] - sample_node->angle_[j];
			if (dis > PI) {
				dis = dis - 2*PI;
			}
			if (dis < -PI) {
				dis = dis + 2*PI;
			}
			sum = sum + dis * dis;
		}
		sum = sqrt(sum);
		if (i == 0) {
			sum_min = sum;
			nearest_node = tree[i]; 
		}
		if (i  > 0) {
			if (sum < sum_min) {
				sum_min = sum;
				nearest_node = tree[i];
			}
		}
	}
	// New config
	double *new_config = new double[numofDOFs];
	// vector from nearest node to sample node;
	double vec_n2s[numofDOFs];
	// vector from nearest node to eps;
	double vec_n2e[numofDOFs];
	double vec[numofDOFs];
	double dis = 0;
	for (int i = 0; i < numofDOFs; i++) {
		vec_n2s[i] = sample_node->angle_[i] - nearest_node->angle_[i];
			if (vec_n2s[i] > PI) {
				vec_n2s[i] = vec_n2s[i] - 2*PI;
			}
			if (vec_n2s[i] < -PI) {
				vec_n2s[i] = vec_n2s[i] + 2*PI;
			}
		dis = dis + vec_n2s[i] * vec_n2s[i];
	}
	dis = sqrt(dis);
	if (dis > eps) {
		// Extend by eps
		double k = eps/dis;
		for (int i = 0; i < numofDOFs; i++) {
			vec_n2e[i] = k * vec_n2s[i];
			new_config[i] = nearest_node->angle_[i] + vec_n2e[i];
		}
		use_eps = true;
	} else {
		// Extend to sample node;
		for (int i = 0; i < numofDOFs; i++) {
			new_config[i] = sample_node->angle_[i];
		}
		use_near = true;
	}
	// check if new node is valid
	if (isObstacleFree(nearest_node->angle_, new_config, numofDOFs, map, x_size, y_size)) {
		node* new_node = new node(numofDOFs);
		for (int i = 0; i < numofDOFs; i++) {
			new_node->angle_[i] = new_config[i];
		}
		new_node->parent_ = nearest_node;
		cost = Distance_between(new_node->angle_, nearest_node->angle_, numofDOFs);
		// cost = MIN(dis, eps);
		// cout << cost << endl;
		new_node->g_ = nearest_node->g_ + cost;
		tree.push_back(new_node);

		if (use_eps) {
			extend_status = 2;
		}
		if (use_near) {
			extend_status = 1;
		}
	} else {
		extend_status = 0;
	}
	return extend_status;
}

int Connect(std::vector<node*> &tree, node* sample_node, int numofDOFs, double* map, int x_size, int y_size, int* n)
{
	int s;
	do {
		s = Extend(tree, sample_node, numofDOFs, map, x_size, y_size);
		n++;
		// cout << s << endl;
	} while(s == 2);

	// cout << "what?" << endl;
	return s;
} 
// update child g value
void update(node_star* node, int numofDOFs) 
{
	cost = Distance_between(node->angle_, node->parent_->angle_, numofDOFs);
	node->g_ = node->parent_->g_ + cost;
	if (node->child_.empty()) {
		return;
	}
	for (int i = 0; i < node->child_.size(); i++) {
		update(node->child_[i], numofDOFs);
	}
}
int Extend_star(std::vector<node_star*> &tree, node_star* sample_node, int numofDOFs, double* map, int x_size, int y_size)
{
	// State for Extend
	bool use_near = false;
	bool use_eps = false;
	double tree_volume;
	int extend_status = 0;//  0 for trapped, 1 for reached, 2 for advanced
	// Get the nearest node
	node_star* nearest_node;
	double dis_s = 0;
	double sum_min = 0;
	double r;
	for (int i = 0; i < tree.size(); i++) {
		double dis = 0;
		double sum = 0;
		for (int j = 0; j < numofDOFs; j++) {
			dis = tree[i]->angle_[j] - sample_node->angle_[j];
			if (dis > PI) {
				dis = dis - 2*PI;
			}
			if (dis < -PI) {
				dis = dis + 2*PI;
			}
			sum = sum + dis * dis;
		}
		sum = sqrt(sum);
		if (i == 0) {
			sum_min = sum;
			nearest_node = tree[i]; 
		}
		if (i  > 0) {
			if (sum < sum_min) {
				sum_min = sum;
				nearest_node = tree[i];
			}
		}
	}
	// New config
	double *new_config = new double[numofDOFs];
	// vector from nearest node to sample node;
	double vec_n2s[numofDOFs];
	// vector from nearest node to eps;
	double vec_n2e[numofDOFs];
	double vec[numofDOFs];
	double dis = 0;
	for (int i = 0; i < numofDOFs; i++) {
		vec_n2s[i] = sample_node->angle_[i] - nearest_node->angle_[i];
			if (vec_n2s[i] > PI) {
				vec_n2s[i] = vec_n2s[i] - 2*PI;
			}
			if (vec_n2s[i] < -PI) {
				vec_n2s[i] = vec_n2s[i] + 2*PI;
			}
		dis = dis + vec_n2s[i] * vec_n2s[i];
	}
	dis = sqrt(dis);
	if (dis > eps) {
		// Extend by eps
		double k = eps/dis;
		for (int i = 0; i < numofDOFs; i++) {
			vec_n2e[i] = k * vec_n2s[i];
			new_config[i] = nearest_node->angle_[i] + vec_n2e[i];
		}
		use_eps = true;
	} else {
		// Extend to sample node;
		for (int i = 0; i < numofDOFs; i++) {
			new_config[i] = sample_node->angle_[i];
		}
		use_near = true;
	}
	// ******************************************************************************
	// check if new node is valid
	if (isObstacleFree(nearest_node->angle_, new_config, numofDOFs, map, x_size, y_size)) {
		node_star* new_node = new node_star(numofDOFs);
		for (int i = 0; i < numofDOFs; i++) {
			new_node->angle_[i] = new_config[i];
		}
		// new_node->parent_ = nearest_node;
		// new_node->g_ = nearest_node->g_ + cost;
		// new_node->parent_->child.push_back(new_node);
		tree.push_back(new_node);
		if (use_eps) {
			extend_status = 2;
		}
		if (use_near) {
			extend_status = 1;
		}
		node_star* min_node = nearest_node;
		
		// rewiring*************************************************
		// Calculate neighborhood radius
		tree_volume = tree.size();
		r =  MIN(pow(gam/del*log(tree_volume)/tree_volume, 1/numofDOFs), eps);	
		// r = PI;
		//************************
		// find the least cost node in the near
		// search through neighbor
		for (int i = 0; i < tree.size(); i++) {
			double d;
			node_star* near_handle = tree[i];
			for (int j = 0; j < numofDOFs; j++) {
				d = d + near_handle->angle_[j] * near_handle->angle_[j];
			}
			d = sqrt(d);
			if (d < r) {
				if (isObstacleFree(near_handle->angle_, new_node->angle_, numofDOFs, map, x_size, y_size)) {
					cost = Distance_between(near_handle->angle_, new_node->angle_,numofDOFs);
					double tmp_g = near_handle->g_ + cost;
					if (tmp_g < min_node->g_) {
						min_node = near_handle;
					}

				}
			}
		}
		// wire new node to least cost node
		new_node->parent_ = min_node;
		cost = Distance_between(new_node->angle_, min_node->angle_, numofDOFs);
		// cout << "cost: " << cost << endl;
		new_node->g_ = min_node->g_ + cost;
		min_node->child_.push_back(new_node);
		node_star* parent_handle;
		// 2nd part in re-wiring
		for(int i = 0; i < tree.size() - 1; i++) {
			node_star* near_handle2 = tree[i];
			double d;
			// cout << "1" << endl;
			// cout << near_handle2 << endl;
			// cout << min_node << endl;
			// cout << endl;
			// cout << tree.back() << endl;
			// cout << tree[tree.size() - 2] << endl;
			// cout << (near_handle2 == min_node) << endl;
			// cout << "*******" << endl;
			if (near_handle2 != min_node) {
				for (int j = 0; j < numofDOFs; j++) {
					d = d + near_handle2->angle_[j] * near_handle2->angle_[j];
				}
				d = sqrt(d);
				// cout << "r: " << r << endl;
				// cout << "d: " << d << endl;
				if (d < r) {
					// cout << "2" << endl;
					if (isObstacleFree(new_node->angle_, near_handle2->angle_, numofDOFs, map, x_size, y_size)) {
						// cout << "3" << endl;
						cost = Distance_between(new_node->angle_, near_handle2->angle_, numofDOFs);
						if (near_handle2->g_ > new_node->g_ + cost) {
							near_handle2->g_ = new_node->g_ + cost;
							parent_handle = near_handle2->parent_;
							// new edge
							near_handle2->parent_ = new_node;
							new_node->child_.push_back(near_handle2);
							// delete previous edge
							for (int j = 0; j < parent_handle->child_.size(); i++) {
								if (parent_handle->child_[j] == near_handle2) {
									parent_handle->child_.erase(parent_handle->child_.begin() + j);
									break;
								}
							}
							// update following g value;
							update(near_handle2, numofDOFs);
							// cout << "odoke" << endl;
						}
					}
				}
			}
		}
	} else {
		extend_status = 0;
	}
	//****************************************************************************
	return extend_status;
}

void Configuarion_generator(int numofDOFs, double* map, int x_size, int y_size)
{
	size_t seed1 = time(NULL);
	std::default_random_engine generator(seed1);
	std::uniform_real_distribution<double> distribution_1(0.0, 2*PI);
	std::uniform_real_distribution<double> distribution_2(0.0, PI);
	double config[numofDOFs];
	do {
		config[0] = distribution_2(generator);
		for (int i = 1; i < numofDOFs; i++) {
			config[i] = distribution_1(generator);
		}
	} while(!IsValidArmConfiguration(config, numofDOFs, map, x_size, y_size));
	cout << "start: " << config[0] << " " << config[1] << " " << config[2] << " " << config[3] << " " << config[4] << endl;
	do {
		config[0] = distribution_2(generator);
		for (int i = 1; i < numofDOFs; i++) {
			config[i] = distribution_1(generator);
		}
	} while(!IsValidArmConfiguration(config, numofDOFs, map, x_size, y_size));
	cout << "goal: " << config[0] << " " << config[1] << " " << config[2] << " "<< config[3] << " " << config[4] << endl;
	return;
}

static void planner(
		   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	//no plan by default
	*plan = NULL;
	*planlength = 0;
    
    //for now just do straight interpolation between start and goal checking for the validity of samples
    double distance = 0;
    int i,j;
    for (j = 0; j < numofDOFs; j++){
        if(distance < fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]))
            distance = fabs(armstart_anglesV_rad[j] - armgoal_anglesV_rad[j]);
    }
    int numofsamples = (int)(distance/(PI/20));
    if(numofsamples < 2){
        printf("the arm is already at the goal\n");
        return;
    }
    *plan = (double**) malloc(numofsamples*sizeof(double*));
    int firstinvalidconf = 1;
    for (i = 0; i < numofsamples; i++){
        (*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
        for(j = 0; j < numofDOFs; j++){
            (*plan)[i][j] = armstart_anglesV_rad[j] + ((double)(i)/(numofsamples-1))*(armgoal_anglesV_rad[j] - armstart_anglesV_rad[j]);
        }
        if(!IsValidArmConfiguration((*plan)[i], numofDOFs, map, x_size, y_size) && firstinvalidconf)
        {
            firstinvalidconf = 1;
            printf("ERROR: Invalid arm configuration!!!\n");
        }
    }    
    *planlength = numofsamples;

    return;
}
static void plannerRRT(
			   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	*plan = NULL;
	*planlength = 0;

	double goal_prob = 0.2;

	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

	size_t seed1 = time(NULL);
	std::default_random_engine generator(seed1);
	std::uniform_real_distribution<double> distribution_1(0.0, 2*PI);
	std::uniform_real_distribution<double> distribution_2(0.0, PI);
	std::uniform_real_distribution<double> distribution_3(0.0, 1.0);

	node* start_node = new node(numofDOFs);
	for (int i = 0; i < numofDOFs; i++) {
		start_node->angle_[i] = armstart_anglesV_rad[i];
	}
	start_node->g_ = 0;
	std::vector<node*> tree;
	tree.push_back(start_node);
	node* sample_node;
	node* new_node;

	double* random_config;
	int goal_flag = 0;
	int n = 0;
	while(true) {
		// Random config
		goal_flag = 0;
		double random_tag = distribution_3(generator);
		if (random_tag < goal_prob) {
			sample_node = new node(numofDOFs);
			for (int i = 0; i < numofDOFs; i ++) {
				sample_node->angle_[i] = armgoal_anglesV_rad[i];
			}
			goal_flag = 1;
		} else {
			random_config = new double[numofDOFs];
			do {
				random_config[0] = distribution_2(generator);
				for (int i = 1; i < numofDOFs; i++) {
					random_config[i] = distribution_1(generator);
				}
			} while(!IsValidArmConfiguration(random_config,numofDOFs, map, x_size, y_size));
			sample_node = new node(numofDOFs);
			for (int i = 0; i < numofDOFs; i++) {
				sample_node->angle_[i] = random_config[i];
			}
		}
		n++;
		// cout << n << endl;
		int extend_status = Extend(tree, sample_node, numofDOFs, map, x_size, y_size);
		if (extend_status == 1 && goal_flag == 1) {
			cout << "Total configuration:  " << n << endl;
			cout << "Total cost: " << tree.back()->g_ << endl;
			break;
		}
		delete sample_node;

	}
	// Find a path;
	node* tmp_node = tree.back();
	int num = 0;
	std::vector<vector<double>> tmp_angle;
	while(tmp_node->parent_ != NULL) {
		std::vector<double> tmp_vec;
		for (int i = 0; i < numofDOFs; i++) {
			tmp_vec.push_back(tmp_node->angle_[i]);
		}
		tmp_angle.push_back(tmp_vec);
		tmp_node = tmp_node->parent_;
		num++;
	}
	int size = tmp_angle.size() + 1;
	*plan = (double**) malloc(size*sizeof(double*));
	// Insert start node
    (*plan)[0] = (double*) malloc(numofDOFs*sizeof(double)); 
    for (int i = 0; i < numofDOFs; i++) {
    	(*plan)[0][i] = armstart_anglesV_rad[i];
    }
	// Insert rest of the node
	for (int i = 0; i < size - 1; i++) {
		(*plan)[i+1] = (double*) malloc(numofDOFs*sizeof(double)); 
		for (int j = 0; j < numofDOFs; j++) {
			(*plan)[i+1][j] = tmp_angle.back()[j];
		}
		tmp_angle.pop_back();
	}
	chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_in_searching = chrono::duration_cast<chrono::duration<double>>(t2-t1);
    printf("time used in searching:  %fs\n",time_in_searching);
	// cout << (*plan)[size-1][0] << " " << (*plan)[size-1][1] << " " << (*plan)[size-1][2] << " " << (*plan)[size-1][3] << " " << (*plan)[size-1][4] << endl;
	*planlength = size;
	// for (int i = 0; i < )
	return;
}

static void plannerRRT_Connect(
		   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

	*plan = NULL;
	*planlength = 0;

	double goal_prob = 0.2;
	size_t seed1 = time(NULL);
	std::default_random_engine generator(seed1);
	std::uniform_real_distribution<double> distribution_1(0.0, 2*PI);
	std::uniform_real_distribution<double> distribution_2(0.0, PI);
	std::uniform_real_distribution<double> distribution_3(0.0, 1.0);

	// Start node
	node* start_node = new node(numofDOFs);
	for (int i = 0; i < numofDOFs; i++) {
		start_node->angle_[i] = armstart_anglesV_rad[i];
	}
	// Goal node
	node* goal_node = new node(numofDOFs);
	for (int i = 0; i < numofDOFs; i++) {
		goal_node->angle_[i] = armgoal_anglesV_rad[i];
	}
	// Start tree
	std::vector<node*> tree_start;
	tree_start.push_back(start_node);
	// Goal tree
	std::vector<node*> tree_goal;
	tree_goal.push_back(goal_node);
	// Swap status flag
	bool extend_start = true;
	bool extend_goal = false;
	node* sample_node;
	node* new_node;

	double* random_config;
	int goal_flag = 0;
	int n = 0;
	while(true) {
		// Random config
		goal_flag = 0;
		double random_tag = distribution_3(generator);
		if (random_tag < goal_prob) {
			sample_node = new node(numofDOFs);
			for (int i = 0; i < numofDOFs; i ++) {
				sample_node->angle_[i] = armgoal_anglesV_rad[i];
			}
			goal_flag = 1;
		} else {
			random_config = new double[numofDOFs];
			do {
				random_config[0] = distribution_2(generator);
				for (int i = 1; i < numofDOFs; i++) {
					random_config[i] = distribution_1(generator);
				}

			} while(!IsValidArmConfiguration(random_config,numofDOFs, map, x_size, y_size));
			// if(n <= 20)
			// cout << n << " " << random_config[0] << " " << random_config[1] << " " << random_config[2] << " " << random_config[3] << " " << random_config[4] << endl;

			sample_node = new node(numofDOFs);
			for (int i = 0; i < numofDOFs; i++) {
				sample_node->angle_[i] = random_config[i];
			}
		}
		n++;
		// Extend start tree;
		if (extend_start == true) {
			int extend_status = Extend(tree_start, sample_node, numofDOFs, map, x_size, y_size);
			n++;
			if (extend_status != 0) {
				new_node = tree_start.back(); 
				int connect_status = Connect(tree_goal, new_node, numofDOFs, map, x_size, y_size, &n);
				if(connect_status == 1) {
					break;
				}
			}
			// Swap
		    extend_start = false;
		    extend_goal = true; 
		}
		if (extend_goal == true) {
			int extend_status = Extend(tree_goal, sample_node, numofDOFs, map, x_size, y_size);
			n++;
			if (extend_status != 0) {
				new_node = tree_goal.back(); 
				int connect_status = Connect(tree_start, new_node, numofDOFs, map, x_size, y_size, &n);
				if(connect_status == 1) {
					break;
				}
			}
			// Swap
			extend_start = true;
			extend_goal = false;
		}
	}

	cout << "connected!" << endl;
	cout << "Total configuration: " << n << endl;
	// Find a path;
	std::vector<vector<double>> tmp_start;
	// THe part in start tree
	node* tmp_node = tree_start.back();
	while(tmp_node->parent_ != NULL) {
		std::vector<double> tmp_sa;
		for (int i = 0; i < numofDOFs; i++) {
			tmp_sa.push_back(tmp_node->angle_[i]);
		}
		// cout << tmp_sa[0] << " " << tmp_sa[1] << " " << tmp_sa[2] << " " << tmp_sa[3] << " " << tmp_sa[4] << " " << endl;
		tmp_start.push_back(tmp_sa);
		tmp_node = tmp_node->parent_;
	}
	// The part in goal tree
	tmp_node = tree_goal.back();
	std::vector<vector<double>> tmp_goal;
	while(tmp_node->parent_ != NULL) {
		std::vector<double> tmp_sa;
		for (int i = 0; i < numofDOFs; i++) {
			tmp_sa.push_back(tmp_node->angle_[i]);
		}
		tmp_goal.push_back(tmp_sa);
		// cout << tmp_sa[0] << " " << tmp_sa[1] << " " << tmp_sa[2] << " " << tmp_sa[3] << " " << tmp_sa[4] << " " << endl;
		tmp_node = tmp_node->parent_;
	}
	int size_start = tmp_start.size();
	int size_goal = tmp_goal.size();
	int size = size_start + size_goal + 2;
	*plan = (double**) malloc(size*sizeof(double*));
	// Insert start node
    (*plan)[0] = (double*) malloc(numofDOFs*sizeof(double)); 
    for (int i = 0; i < numofDOFs; i++) {
    	(*plan)[0][i] = armstart_anglesV_rad[i];
    }
    for (int i = 0; i < size_start; i++) {
    	(*plan)[i + 1] = (double*) malloc(numofDOFs*sizeof(double)); 
		for(int j = 0; j < numofDOFs; j++) {
			(*plan)[i + 1][j] = tmp_start.back()[j];
		}
		tmp_start.pop_back();
    }
    for (int i = 0; i < size_goal; i++) {
    	int m = size_start + i;
		(*plan)[m + 1] = (double*) malloc(numofDOFs*sizeof(double)); 
		for(int j = 0; j < numofDOFs; j++) {
			(*plan)[m + 1][j] = tmp_goal[i][j];
		}
    }
	(*plan)[size - 1] = (double*) malloc(numofDOFs*sizeof(double)); 
	for (int i = 0; i < numofDOFs; i ++) {
		(*plan)[size - 1][i] = armgoal_anglesV_rad[i];
	}

	chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_in_searching = chrono::duration_cast<chrono::duration<double>>(t2-t1);
    printf("time used in searching:  %fs\n",time_in_searching);
	*planlength = size;


	// for (int i = 0; i < 20; i++) {
 //    	cout << i << endl;
 //    	Configuarion_generator(numofDOFs, map, x_size,y_size);
 //    }

	return;
}

static void plannerRRT_Star(
		   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();

	*plan = NULL;
	*planlength = 0;
	int n = 0;

	int goal_flag;
	double goal_prob = 0.2;
	size_t seed1 = time(NULL);
	std::default_random_engine generator(seed1);
	std::uniform_real_distribution<double> distribution_1(0.0, 2*PI);
	std::uniform_real_distribution<double> distribution_2(0.0, PI);
	std::uniform_real_distribution<double> distribution_3(0.0, 1.0);
	vector<node_star*> V;
	// start goal
	node_star* node_start = new node_star(numofDOFs);
	for (int i = 0; i < numofDOFs; i++) {
		node_start->angle_[i] = armstart_anglesV_rad[i];
	}
	node_start->g_ = 0;
	V.push_back(node_start);

	node_star* sample_node;

	double* random_config;
	// start the loop
	while(true) {
		// Random config
		goal_flag = 0;
		double random_tag = distribution_3(generator);
		if (random_tag < goal_prob) {
			sample_node = new node_star(numofDOFs);
			for (int i = 0; i < numofDOFs; i ++) {
				sample_node->angle_[i] = armgoal_anglesV_rad[i];
			}
			goal_flag = 1;
		} else {
			random_config = new double[numofDOFs];
			do {
				random_config[0] = distribution_2(generator);
				for (int i = 1; i < numofDOFs; i++) {
					random_config[i] = distribution_1(generator);
				}

			} while(!IsValidArmConfiguration(random_config,numofDOFs, map, x_size, y_size));
			sample_node = new node_star(numofDOFs);
			for (int i = 0; i < numofDOFs; i++) {
				sample_node->angle_[i] = random_config[i];
			}
		}
		n++;
		int extend_flag = Extend_star(V, sample_node, numofDOFs, map, x_size, y_size);
		// if (n == 1000000)
		// break;
		// cout << extend_flag << endl;
		if (goal_flag == 1 && extend_flag == 1) {
			cout << "Path found!" << endl;
			cout << "Total configuration: " << n << endl;
	        cout << "Total cost: " << V.back()->g_ << endl;
			break;
		}
		// cout << "check !" << endl;
	}
	// Find a path;
	node_star* tmp_node = V.back();
	int num = 0;
	std::vector<vector<double>> tmp_angle;
	while(tmp_node->parent_ != NULL) {
		std::vector<double> tmp_vec;
		for (int i = 0; i < numofDOFs; i++) {
			tmp_vec.push_back(tmp_node->angle_[i]);
		}
		tmp_angle.push_back(tmp_vec);
		tmp_node = tmp_node->parent_;
		num++;
	}
	int size = tmp_angle.size() + 1;
	*plan = (double**) malloc(size*sizeof(double*));
	// Insert start node
    (*plan)[0] = (double*) malloc(numofDOFs*sizeof(double)); 
    for (int i = 0; i < numofDOFs; i++) {
    	(*plan)[0][i] = armstart_anglesV_rad[i];
    }
	// Insert rest of the node
	for (int i = 0; i < size - 1; i++) {
		(*plan)[i+1] = (double*) malloc(numofDOFs*sizeof(double)); 
		for (int j = 0; j < numofDOFs; j++) {
			(*plan)[i+1][j] = tmp_angle.back()[j];
		}
		tmp_angle.pop_back();
	}
	chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_in_searching = chrono::duration_cast<chrono::duration<double>>(t2-t1);
    printf("time used in searching:  %fs\n",time_in_searching);

	// cout << (*plan)[size-1][0] << " " << (*plan)[size-1][1] << " " << (*plan)[size-1][2] << " " << (*plan)[size-1][3] << " " << (*plan)[size-1][4] << endl;
	*planlength = size;
	return;
}
static void plannerPRM(
			   double*	map,
		   int x_size,
 		   int y_size,
           double* armstart_anglesV_rad,
           double* armgoal_anglesV_rad,
	   int numofDOFs,
	   double*** plan,
	   int* planlength)
{
	*plan = NULL;
	*planlength = 0;
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();


	int comp_index = 0;
	int N = 40000;
	double radius_neighbor = 0.2*PI;
	size_t seed1 = time(NULL);
	std::default_random_engine generator(seed1);
	std::uniform_real_distribution<double> distribution_1(0.0, 2*PI);
	std::uniform_real_distribution<double> distribution_2(0.0, PI);

	std::vector<node_prm*> G;
	node_prm* sample_node;
	double* random_config;
	int num_of_node = 0;
	int edge_num = 0;
	int n = 0;
	// Build PRM
	while (num_of_node < N) {
		comp_index++;
		num_of_node++;
		// cout << "1" << endl;
		random_config = new double[numofDOFs];
		// Random configuration
		do {
			random_config[0] = distribution_2(generator);
			for (int i = 1; i < numofDOFs; i++) {
				random_config[i] = distribution_1(generator);
			}

		} while(!IsValidArmConfiguration(random_config,numofDOFs, map, x_size, y_size));
		// cout << "2" << endl;
		sample_node = new node_prm(numofDOFs, comp_index);
		for (int i = 0; i < numofDOFs; i++) {
			sample_node->angle_[i] = random_config[i];
		}
		G.push_back(sample_node);
		// Find the neighborhood of sampled node;
		std::vector<node_prm*> neighborhood;
		for (int i = 0; i < G.size() - 1; i++) {
			if (Distance_between(G[i]->angle_, sample_node->angle_, numofDOFs) < radius_neighbor) {
				neighborhood.push_back(G[i]);
			}
		}		
		// Connect newly sampled node to component(neighbor)
		for (int i = 0; i < neighborhood.size(); i++) {
			// cout << "5" << endl;
			if (neighborhood[i]->comp_index_ != sample_node->comp_index_ && isObstacleFree(neighborhood[i]->angle_, sample_node->angle_, numofDOFs, map, x_size,y_size) == 1) {
				int previous_index = MIN(neighborhood[i]->comp_index_, sample_node->comp_index_);
				int later_index = MAX(neighborhood[i]->comp_index_, sample_node->comp_index_);
				// cout << "6" << endl;
				for (int j = 0; j < G.size(); j++) {
					if (G[j]->comp_index_ == later_index) 
						G[j]->comp_index_ = previous_index;
				}
				edge_num ++;
				neighborhood[i]->connected.push_back(sample_node);
				sample_node->connected.push_back(neighborhood[i]);
			}
		}
	}
	chrono::steady_clock::time_point t3 = chrono::steady_clock::now();
    chrono::duration<double> time_in_construction = chrono::duration_cast<chrono::duration<double>>(t3-t1);
    printf("construction time:  %fs\n",time_in_construction);
	// cout << " Total edge: " << edge_num << endl;
	cout << "num of config: " << num_of_node << endl;
	// cout << "n: " << n << endl;
	// Generate start node and goal node
	comp_index++;
	node_prm* node_start = new node_prm(numofDOFs, comp_index);
	comp_index++;
	node_prm* node_goal = new node_prm(numofDOFs, comp_index);
	for (int i = 0; i < numofDOFs; i++) {
		node_start->angle_[i] = armstart_anglesV_rad[i];
		node_goal->angle_[i] = armgoal_anglesV_rad[i];
	}

	// Find the reachable nodes for start and goal node;
	std::vector<node_prm*> reachable_start;
	std::vector<node_prm*> reachable_goal;
	for (int i = 0; i < G.size(); i++) {
		if (Distance_between(G[i]->angle_, node_start->angle_, numofDOFs) < radius_neighbor
			&& isObstacleFree(G[i]->angle_, node_start->angle_, numofDOFs, map, x_size, y_size)) {
			reachable_start.push_back(G[i]);
		}
		if (Distance_between(G[i]->angle_, node_goal->angle_, numofDOFs) < radius_neighbor
			&& isObstacleFree(G[i]->angle_, node_goal->angle_, numofDOFs, map, x_size, y_size)) {
			reachable_goal.push_back(G[i]);
		}
	}

	// Find the component that both goal and start can connect
	int component_index = 0;
	int reachable_start_index = -1;
	int reachable_goal_index = -1;
	// Search through the start reachable set
	for (int i = 0; i < reachable_start.size(); i++) {
		// cout << "3" << endl;
		// If start node can be connected to some component
		if (isObstacleFree(reachable_start[i]->angle_, node_start->angle_, numofDOFs, map, x_size, y_size)) {
			// Search through the goal reachable set
			// cout << "4" << endl;
			for (int j = 0; j < reachable_goal.size(); j++) {
				// cout << "5" << endl;
				// If goal node can be connected to some component
				if (isObstacleFree(reachable_goal[j]->angle_, node_goal->angle_, numofDOFs, map, x_size, y_size)) {
					// Check if they are same component
					if (reachable_start[i]->comp_index_ == reachable_goal[j]->comp_index_) {
						component_index = reachable_start[i] ->comp_index_;
						reachable_start_index = i;
						reachable_goal_index = j;
						break;
					}
				}
			}
			if (component_index != 0)
				break;
		}
	}

	// Break if no path is found;
	if(component_index == 0) {
		cout << "No path found." << endl;
		return;
	}

	
	// Add start and goal node into the RPM
	reachable_start[reachable_start_index]->connected.push_back(node_start);
	node_start->connected.push_back(reachable_start[reachable_start_index]);
	reachable_goal[reachable_goal_index]->connected.push_back(node_goal);
	node_goal->connected.push_back(reachable_goal[reachable_goal_index]);

	// A-star search from start to goal
	// build and initialize Open list
	vector<node_prm*> open;
	open.push_back(node_start);
	make_heap(open.begin(), open.end(), compare_func());
	// build the closed list;
	unordered_map<int, node_prm*> closed;
	while(!open.empty()) {
		node_prm* current_node = open.front();
		int current_index = current_node->node_index_;
		pop_heap(open.begin(), open.end(), compare_func());
		open.pop_back();
		closed.insert(pair<int, node_prm*>(current_index, current_node));

		for (int i = 0; i < current_node->connected.size(); i++) {
			// check if is in clsoed 
			node_prm* successor = current_node->connected[i];
			std::unordered_map<int, node_prm*>::const_iterator it = closed.find(successor->node_index_);
			if (it == closed.end()) {
				// cout << "douo" << endl;
				// check if is in open
				bool found_in_open = false;
				for (int j = 0; j < open.size(); j++) {
					// if found in open
					if (open[j]->node_index_ == successor->node_index_) {
						found_in_open = true;
						// check the g value
						double cost = Distance_between(open[j]->angle_, current_node->angle_, numofDOFs);
						if (open[j]->g_ > current_node->g_ + cost) {
							open[j]->g_ = current_node->g_ + cost;
							open[j]->parent_ = current_node;
						}
						break;
					}
				}
				if (!found_in_open) {
					double cost = Distance_between(current_node->angle_, successor->angle_, numofDOFs);
					successor->g_ = current_node->g_ + cost;
					successor->parent_ = current_node;
					open.push_back(successor);
					push_heap(open.begin(), open.end(), compare_func());
				}		
			}
		}
		if(current_index == node_goal->node_index_) {
			cout << "Path found" << endl;
			break;
		}
	}

	// back track
	std::unordered_map<int, node_prm*>::const_iterator it = closed.find(node_goal->node_index_);
	std::vector<vector<double>> angle_vec;
	
	cout << it->second->angle_[1] << endl;
	if(it != closed.end()){
		// cout << "enen" << endl;
		node_prm* handle = it->second;
		do {
			std::vector<double> angle_tmp;
			for (int i = 0; i < numofDOFs; i++) {
				angle_tmp.push_back(handle->angle_[i]);
			}
			angle_vec.push_back(angle_tmp);
			handle = handle->parent_;
		} while (handle != NULL);
	}
	int size = angle_vec.size();
	*plan = (double**) malloc(size*sizeof(double*));

	for (int i = 0; i < size; i++) {
		(*plan)[i] = (double*) malloc(numofDOFs*sizeof(double)); 
		for (int j = 0; j < numofDOFs; j++) {
			(*plan)[i][j] = angle_vec.back()[j];
		}
		angle_vec.pop_back();
	}
	chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
    chrono::duration<double> time_in_searching = chrono::duration_cast<chrono::duration<double>>(t2-t1);
    printf("time used in searching:  %fs\n",time_in_searching);
    cout << node_goal->g_ << endl;
	*planlength = size;
	return;
}
//prhs contains input parameters (3): 
//1st is matrix with all the obstacles
//2nd is a row vector of start angles for the arm 
//3nd is a row vector of goal angles for the arm 
//plhs should contain output parameters (2): 
//1st is a 2D matrix plan when each plan[i][j] is the value of jth angle at the ith step of the plan
//(there are D DoF of the arm (that is, D angles). So, j can take values from 0 to D-1
//2nd is planlength (int)
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[])
     
{ 
    
    /* Check for proper number of arguments */    
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidNumInputs",
                "Four input arguments required."); 
    } else if (nlhs != 2) {
	    mexErrMsgIdAndTxt( "MATLAB:planner:maxlhs",
                "One output argument required."); 
    } 
        
    /* get the dimensions of the map and the map matrix itself*/     
    int x_size = (int) mxGetM(MAP_IN);
    int y_size = (int) mxGetN(MAP_IN);
    double* map = mxGetPr(MAP_IN);
    
    /* get the start and goal angles*/     
    int numofDOFs = (int) (MAX(mxGetM(ARMSTART_IN), mxGetN(ARMSTART_IN)));
    if(numofDOFs <= 1){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "it should be at least 2");         
    }
    double* armstart_anglesV_rad = mxGetPr(ARMSTART_IN);
    if (numofDOFs != MAX(mxGetM(ARMGOAL_IN), mxGetN(ARMGOAL_IN))){
        	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidnumofdofs",
                "numofDOFs in startangles is different from goalangles");         
    }
    double* armgoal_anglesV_rad = mxGetPr(ARMGOAL_IN);
 
    //get the planner id
    int planner_id = (int)*mxGetPr(PLANNER_ID_IN);
    if(planner_id < 0 || planner_id > 3){
	    mexErrMsgIdAndTxt( "MATLAB:planner:invalidplanner_id",
                "planner id should be between 0 and 3 inclusive");         
    }
    
    //call the planner
    double** plan = NULL;
    int planlength = 0;
    
    //you can may be call the corresponding planner function here
    if (planner_id == RRT)
    {
       plannerRRT(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }

    if (planner_id == RRTCONNECT)
    {
       plannerRRT_Connect(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }
    if (planner_id == RRTSTAR) 
    {
	   plannerRRT_Star(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }
    if (planner_id == PRM) 
    {
    	plannerPRM(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength);
    }
    


    //dummy planner which only computes interpolated path

    // planner(map,x_size,y_size, armstart_anglesV_rad, armgoal_anglesV_rad, numofDOFs, &plan, &planlength); 
    
    printf("planner returned plan of length=%d\n", planlength); 
    
    /* Create return values */
    if(planlength > 0)
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)planlength, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);        
        //copy the values
        int i,j;
        for(i = 0; i < planlength; i++)
        {
            for (j = 0; j < numofDOFs; j++)
            {
                plan_out[j*planlength + i] = plan[i][j];
            }
        }
    }
    else
    {
        PLAN_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)numofDOFs, mxDOUBLE_CLASS, mxREAL); 
        double* plan_out = mxGetPr(PLAN_OUT);
        //copy the values
        int j;
        for(j = 0; j < numofDOFs; j++)
        {
                plan_out[j] = armstart_anglesV_rad[j];
        }     
    }
    PLANLENGTH_OUT = mxCreateNumericMatrix( (mwSize)1, (mwSize)1, mxINT8_CLASS, mxREAL); 
    int* planlength_out = (int*) mxGetPr(PLANLENGTH_OUT);
    *planlength_out = planlength;

    
    return;
    
}





