/*
 * A skeleton program for a branch-and-bound algorithm for the TSP
 * */

#ifndef DEBUG
#define DEBUG 0
#endif

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>
#include <stack>
#include <cstdlib>
#include <map>
#include <cmath>
#include <climits>
#include <random>
#include <exception>
#include <time.h>
#include <sys/time.h>
#include <chrono>

using namespace std;
using namespace std::chrono;


typedef size_t Vertex;
typedef pair <Vertex, Vertex> Edge;
typedef map <Vertex, vector <Edge> > Incidence_List;
typedef vector <Edge>::iterator edge_vector_iter;
typedef map <Edge, double> Weight_Map;

struct timeval tv_global;

ostream& operator <<(ostream &os, const Edge &e){
        os << "(" << e.first << "," << e.second << ")";
        return os;
}

class Graph {
public:

set <Vertex> vertices;
set <Edge> edges;
Incidence_List incList;
Weight_Map weightMap;
bool isDirected;

// Let us have a simple default constructor
Graph() :
        vertices(),
        edges(),
        incList(),
        weightMap(),
        isDirected(false)
{
        // pass
}


size_t numOfVerts(){
        return this->vertices.size();
}

size_t numOfEdges(){
        return this->edges.size();
}

/*
 * a helping function to tell us the degree of a vertex
 */
size_t degree(const Vertex &v){
        try{
                size_t deg = this->incList[v].size();
                return deg;
        }
        catch (exception &e) {
                cerr << "Something happened in getting the degree of "
                     << v <<". Maybe it is not in the graph?" << endl;
                return -1;
        }
}

/*
 * A function to add a weighted edge.
 * The default weight of an edge is set to be 0.
 * The function will add a single directed edge
 * if you want your graph to behave as undirected,
 * do not forget to add the same edge in the opposite direction
 */
void addDirectedEdge(const Edge &e, double w = 0){
        // What if the edge already exists?
        if (this->edges.find(e) != this->edges.end()) {
//			cerr << "The egde " << e << " already exists in the graph!" << endl;
                return;
        }

        Vertex u, v;
        u = e.first;
        v = e.second;

        /*
         * Don't forget to add the vertices if they do not already exist!
         */
        this->vertices.insert(u);
        this->vertices.insert(v);

        /*
         * Same goes for the edge e
         */
        this->edges.insert(e);

        /*
         * Let us map the weight of edge e
         */
        this->weightMap[e] = w;

        /*
         * Finally, we put the edge e into the incidence list of
         * vertex u
         */
        try{
                this->incList[u].push_back(e);
        }
        catch (exception &ex) {
                cerr << "Error, something happened when trying to "
                     << "append edge " << e << endl;
        }

        return;
}

/* Let us overload the addDirectedEdge method
 * to take two vertices instead an edge
 */
void addDirectedEdge(const Vertex &u, const Vertex &v, double w = 0){
        // we will just turn the vertex pair (u,v) into an edge
        // and then use the function which takes an edge as an argument
        Edge e = make_pair(u, v);
        this->addDirectedEdge(e, w);
        return;
}

/*
 * It might be useful to remove edges from graphs
 * Be careful to remove the edge directed oppositely if the graph is not directed
 */
void removeDirectedEdge(const Edge &e){
        // Let us check if the edge e is registered in the
        // graph's set of edges
        if (this->edges.find(e) == this->edges.end()) {
                cerr << "Strange, the edge " << e
                     <<" is not registered in the graph!" << endl;
        } // if

        // erasing from sets is not that difficult
        try{
                this->edges.erase(e);
        }
        catch (exception &ex) {
                cerr << "Error erasing edge " << e << endl;
        }

        // neither is from maps
        this->weightMap.erase(e);

        // erasing from vectors is a different story
        Vertex u;
        u = e.first;

        // Erase e from u's incidence list
        auto u_begin = this->incList[u].begin();
        auto u_end = this->incList[u].end();
        this->incList[u].erase(remove(u_begin, u_end, e), u_end);

//		/*
//		 * Now, if the graph is not directed, and the degree of u becomes 0,
//		 * then we will also remove u from the list of vertices.
//		 * This may not be good if the graph is directed, because there may be
//		 * some edge pointing to it!
//		*/
//		if (not this->isDirected){
//			if (this->degree(u) == 0){
//				this->vertices.erase(u);
//			}
//		}
        return;
}     // removeDirectedEdge

/* Let's overload the edge removal function
 * to work with two given vertices
 */
void removeDirectedEdge(const Vertex& u, const Vertex& v){
        Edge e = make_pair(u, v);
        this->removeDirectedEdge(e);
}


/*
 * A helpful function which will print the graph's incidence list for us
 */
void print_incList(){
        for (auto v_iter = this->vertices.begin();
             v_iter != this->vertices.end();
             v_iter++) {
                Vertex v = *v_iter;
                cout << v << ": ";
                for (auto e_iter = this->incList[v].begin();
                     e_iter!= this->incList[v].end();
                     e_iter++) {
                        Edge e = *e_iter;
                        cout << "[" << e << "," << this->weightMap[e] << "], ";
                } // for all incident edges to v
                cout << endl;
        } // for all vertices v
}     // void print_incList()

/*
 * A DFS (depth-first-search) function
 * takes a starting vertex and
 * stores the predecessor of each vertex in a given map
 */
void DFS(Vertex s, map<Vertex, Vertex> &preds){
        stack <Vertex> d_stack;
        map<Vertex, bool> discovered;

        for (auto v_iter = this->vertices.begin();
             v_iter != this->vertices.end();
             v_iter++) {
                Vertex v = *v_iter;
                discovered[v] = false;
                preds[v] = v; // Initialize each vertex as its own predecessor
        }

        d_stack.push(s);

        while (not d_stack.empty()) {
                Vertex v = d_stack.top();
                d_stack.pop();
                discovered[v] = true;
                for (auto e_iter = this->incList[v].begin();
                     e_iter != this->incList[v].end();
                     e_iter++) {
                        Edge e = *e_iter;
                        Vertex u = e.second;
                        if (not discovered[u]) {
                                discovered[u] = true;
                                preds[u] = v;
                                d_stack.push(u);
                        } // if not discovered[u]
                } // for incident edges to v
        } // while d_stack is not empty
}     // function DFS

/*
 * We would like to know if a graph is a cycle
 */
bool is_cycle(){
        // Let's check if the graph has vertices
        if (this->numOfVerts() == 0) {
                cerr << "The graph is empty!" <<endl;
                return false;
        }

        /* We will use dfs to detect if a graph is a cycle
         * !! Be careful:
         * We want to check if the complete graph is a cycle,
         * not just if it contains a cycle or not.
         * For this, we will use DFS
         */
        // Let us initialize a predecessor map
        map <Vertex, Vertex> preds;

        stack <Vertex> d_stack;
        map<Vertex, bool> processed;

        for (auto v_iter = this->vertices.begin();
             v_iter != this->vertices.end();
             v_iter++) {
                Vertex v = *v_iter;
                processed[v] = false;
                preds[v] = v; // Initialize each vertex as its own predecessor
        }

        // Then, take any vertex of the graph
        Vertex s = *(this->vertices.begin());
        d_stack.push(s);

        while (not d_stack.empty()) {
                Vertex v = d_stack.top();
                d_stack.pop();
                if (not processed[v]) {
                        processed[v] = true;
                        for (auto e_iter = this->incList[v].begin();
                             e_iter != this->incList[v].end();
                             e_iter++) {
                                Edge e = *e_iter;
                                Vertex u = e.second;
                                if (not processed[u]) {
                                        // discovered[u] = true;
                                        preds[u] = v;
                                        d_stack.push(u);
                                } // if not discovered[u]
                                else{
                                        if (preds[v] != u) {
                                                if (u != s) {
                                                        /* This means that we have found a backward edge
                                                         * that does not go to the initial vertex s
                                                         */
                                                        return false;
                                                } // if (preds[v] != u)
                                                else break; // stop the dfs traversal
                                        }
                                } // else (u has been discovered before)
                        } // for incident edges to v
                } // if (not discovered[v])
        } // while d_stack is not empty

        // Let us check if all vertices were processed by the DFS
        for (auto v_iter = this->vertices.begin();
             v_iter != this->vertices.end();
             v_iter++) {
                Vertex u = *v_iter;
                if (not processed[u]) {
                        return false;
                } // if
        } // for all vertices

        /* As final check, we check if the degree
         * of each vertex is equal to feasible_degree!
         */
        size_t feasible_degree = 2;
        if (this->isDirected) feasible_degree = 1;
        for (auto v_iter = this->vertices.begin();
             v_iter != this->vertices.end();
             v_iter++) {
                Vertex u = *v_iter;
                if (this->degree(u) != feasible_degree) {
                        return false;
                }     //if
        } // for all vertices

        /* So far nothing strange has happened.
         * The graph must be a cycle.
         */
        return true;
}     // is_cycle()

double total_weight(){
        double sum = 0;
        for (auto e_iter = this->edges.begin();
             e_iter != this->edges.end();
             e_iter++) {
                Edge e = *e_iter;
                double add = this->weightMap[e];
                if (not this->isDirected) add /= 2.0;
                sum += add;
        } // for all edges
        return sum;
}
};

double timedif(const struct timeval& tv_start, const struct timeval& tv_end){
        unsigned long seconds, useconds;
        seconds = tv_end.tv_sec - tv_start.tv_sec - 1;
        useconds = 1000000 + tv_end.tv_usec - tv_start.tv_usec;
        return (double)seconds + (double)useconds/1000000.0;
}

double forcedTSP(   Graph &G,
                    Graph &F
                    /* any parameters related to bounding? */
                    )
{
        // What if G and F are different types directed/undirected?
        if (G.isDirected != F.isDirected) {
                cerr << "Why are G and F different?" << endl;
                return INFINITY;
        }


        static unsigned long long branches = 0;
        static unsigned long long leaves = 0;
        static double incumbent = INFINITY;
        //You can use the "incumbent" variable for bounding operations!

        branches++; // Let's just keep track how many times we've branched
        if (DEBUG) {
                // Print some messages to track the progress
                if (not (branches % 50000)) {
                        struct timeval tv_local;
                        gettimeofday(&tv_local, NULL);
                        cout << "Branches: " << branches << endl;
                        cout << "Elapsed time: " << timedif(tv_global, tv_local)
                             << " seconds" << endl;
                }
        }

        // Let us check if G is empty
        if (G.numOfEdges() == 0) {
                leaves++; // We've reached a leaf in the solution tree!

                // Then, let's see if F is a Hamiltonian cycle in G
                if (F.numOfVerts() == G.numOfVerts() && F.is_cycle()) {
                        double tw = F.total_weight();
                        if (incumbent > tw) incumbent = tw;
                        // we update the invumbent solution
                        return tw;
                }
                else return INFINITY; // return INFINITY if F is not a solution

                if (DEBUG) {
                        // Print some messages to track the progress
                        if (not (leaves % 10000)) {
                                struct timeval tv_local;
                                gettimeofday(&tv_local, NULL);
                                cout << "Leaves: " << leaves << endl;
                                cout << "Elapsed time: " << timedif(tv_global, tv_local)
                                     << " seconds" << endl;
                                cout << "Current value: " << incumbent << endl;
                        } // 10000 leaves
                } //if (DEBUG)
        } // if (G.numOfEdges() == 0),
          // that is, this is a leaf of the search tree
        else{

                /********************************************************
                *
                *       Add your own
                *       REDUCTION and BOUNDING operations
                *
                ********************************************************/


                // just pick an edge (is picking any edge good enough?)
                Edge e = *G.edges.begin();

                // prepare 2 copies of F
                Graph F_force(F);
                Graph F_delete(F);

                // A single copy of G is enough
                Graph Gfd(G);

                Gfd.removeDirectedEdge(e);
                F_force.addDirectedEdge(e, G.weightMap[e]);

                if (not Gfd.isDirected) {
                        Edge ee = make_pair(e.second, e.first);
                        Gfd.removeDirectedEdge(ee);
                        F_force.addDirectedEdge(ee, G.weightMap[ee]);
                }

                // call the Force branch recursively
                double ans_force;
                ans_force = forcedTSP(Gfd, F_force);

                // call the Delete branch recursively
                double ans_delete;
                ans_delete = forcedTSP(Gfd, F_delete);

                if (ans_force < ans_delete) {
                        F = Graph(F_force);
                        return ans_force;
                }
                else{
                        F = Graph(F_delete);
                        return ans_delete;
                }
        } //else (G is not empty!)
}

double uniform(double min, double max){
        return min + (max-min)*rand()/(double(RAND_MAX));
}

/* A helper function to help us get rid of whitespace
 * surrounding a string
 */
inline std::string trim(const std::string &s)
{
        auto wsfront=std::find_if_not(s.begin(),s.end(),[](int c){
                return std::isspace(c);
        });
        return std::string(wsfront,std::find_if_not(s.rbegin(),std::string::const_reverse_iterator(wsfront),[](int c){
                return std::isspace(c);
        }).base());
}

/* A function to read a tsplib file from an input stream,
 * default is cin.
 * return value should be:
 * 1 - all ok
 * 0 - something happened, graph not fully read
 */
int read_tsplib(Graph &g, istream& is = cin){
        string estr, strpre, strpost;
        size_t n = 0; // the # vertices of the graph
        char colon = ':';
        const string dimension = string("DIMENSION");
        const string weight_format = string("EDGE_WEIGHT_FORMAT");
        const string weight_type = string("EDGE_WEIGHT_TYPE");
        const string good_weight_type = string("EXPLICIT");
        const string good_weight_format = string("FULL_MATRIX");
        const string edge_weights_follow = string("EDGE_WEIGHT_SECTION");
        const string str_eof = string("EOF");

        while (true) {
                try{
                        getline(is, estr);
                        stringstream strm(estr);
                        getline(strm, strpre, colon);
                        getline(strm, strpost);
                        strpre = trim(strpre); // remove leading and trailing whitespace
                        strpost = trim(strpost);
                }
                catch (exception &e) {
                        cerr << "That was exceptional!" << endl;
                }


                if (not strpre.compare(dimension)) {
                        // The function string.compare() returns 0
                        // if the two strings are equal
                        n = stoi(strpost);
                }
                else if (not strpre.compare(weight_type)) {
//				getline(is, strpost);
//				strpost = trim(strpost);
                        if (strpost.compare(good_weight_type)) {
                                // The function string.compare() returns 0
                                // if the two strings are equal
                                cout << "Please use EXPLICIT edge weights in the "
                                     << "tsplib file." << endl;
                                return 0;
                        }
                }
                else if (not strpre.compare(weight_format)) {
//			getline(is, strpost);
//			strpost = trim(strpost);
                        if (strpost.compare(good_weight_format)) {
                                // The function string.compare() returns 0
                                // if the two strings are equal
                                cout << "Please use FULL_MATRIX format in the "
                                     << "tsplib file." << endl;
                                return 0;
                        }
                }
                else if (not strpre.compare(edge_weights_follow)) {
                        size_t i, j;
                        for (i = 0; i < n; i++)
                                for (j = 0; j < n; j++) {
                                        double weight;
                                        is >> weight;
                                        if (i == j) continue;
                                        Vertex u = i;
                                        Vertex v = j;
                                        g.addDirectedEdge(u, v, weight);
                                }
                        return 1;
                }
                else if (not strpre.compare("str_eof")) {
                        return 0;
                        // reached the end...
                }
                else {
//			getline(is, strpost);
                }
        }
}


int main(int argc, char **argv){
        Graph g;
        g.isDirected = false;
        read_tsplib(g);
        // cout << "The input graph is:" << endl;
        // g.print_incList();
        // cout << endl;

        Graph ff;
        // cout << "Starting calculation..." << endl;
        double solVal, solTime;
        timeval tv_solutionTime;

        gettimeofday(&tv_global, NULL);
        steady_clock::time_point start = steady_clock::now();

        solVal = forcedTSP(g, ff);

        gettimeofday(&tv_solutionTime, NULL);

        steady_clock::time_point end = steady_clock::now();
        steady_clock::duration d = end - start;

        // solTime = timedif(tv_global, tv_solutionTime);
        solTime = duration_cast<milliseconds>(d).count();

        // double LIMINF = 1.84467e+10;
        // if (solTime > LIMINF){
        //   int n = 10000;
        //   gettimeofday(&tv_global, NULL);
        //
        //   for(int i = 0; i < n; i++){
        //     Graph ff;
        //     double lbound = 0;
        //     solVal = forcedTSP(g, ff, lbound);
        //   }
        //   gettimeofday(&tv_solutionTime, NULL);
        //   solTime =  timedif(tv_global, tv_solutionTime) / float(n);
        // }
        // cout << "FINISHED!" << endl;
        // cout << "\tElapsed Time: " << solTime << " seconds." << endl;
        // cout << "Optimal Route Length: " << solVal << endl;
        // cout << "The route:" << endl;
        // ff.print_incList();
        std::cout << solTime << endl;
}
