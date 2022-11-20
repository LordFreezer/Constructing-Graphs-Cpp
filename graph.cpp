#include <iostream>
#include <vector>
#include <queue>

using namespace std;
class Edge
// class to keep track of things used often 
{
public:
    
    // variable for edge count
    int edgeNum;
    
    // variable to keep track of isolated nodes that have been called twice
    int doubleCount; 
    
    // queue to store the isolated vertices
    queue<int> isolatedNodes;
    
    // simple char array that lets me use letters instear of numbers
    char A[14] = { 'A','B','C','D','E','F','G','H','I',
              'J','K','L','M','N' };

    // Create a vector to store a vertex and its connected vertices.
    // This can be used with both a directed graph and undirected graph.
    vector<int> adj[7];

    // Creates a vector that reverses the direction that the nodes connects to in
    // a directed graph. If the user selects a undirected graph, then adj will do
    // the same thing as rev + adj as it connects nodes both ways.
    vector<int> trans[7];

    char type;
} e;

void adjacentVertex(int** matrix, int ** transposed, int V) 
// looks for vertices in a matrix that are mutually adjacent
{
    
    for (int i = 0; i < V - 1; i++) 
    {
        for (int j = 0; j < V - 1; j++) 
        {
            if ((matrix[i][j] == transposed[i][j])&&(matrix[i][j] == 1))
                cout << e.A[i] << " is mutually connected to " << e.A[j] << endl;
        }
    }
}

int min(int L, int R)
{
    // base case
    if (L == R)
        return L;
    // recursive call
    if (L > R)
        return min(L - 1, R);
    // last resort case
    else
        return L;
}

int max(int L, int R)
{
    // base case
    if (L == R)
        return L;
    // recursive call
    if (L < R)
        return max(L + 1, R);
    // last resort case
    else
        return L;
}

int fun(int A[], int L, int R, bool condition)
// recursivley divides an array to find its extrema
{
    // base case
    if (L == R)
        return (A[L]);
    int M = (L + R) / 2;
    // recursive tail call
    switch (condition)
    {
    case true:
        return max(fun(A, L, M, condition), fun(A, M + 1, R, condition));
    case false:
        return min(fun(A, L, M, condition), fun(A, M + 1, R, condition));
    }
}

void Degree(int** matrix, int V, char type) 
{
    // This lets me reuse code
    if(type == 'o')
        cout << "OUT DEGREE" << endl;
    else
        cout << "IN DEGREE" << endl;

    int count; // This is actually the number of vertices connected to the each vertex 
    for (int i = 0; i < V - 1; i++) 
    {
        count = 0;
        // generates number of vertices
        for (int j = 0; j < V - 1; j++) 
        {
            if (matrix[i][j] > -1) 
            {
                count++;
            }          
        }
        // index used for storing vertices
        int blob = 0;
        // dynamically allocated array that stores vertices.
        // size is based on which vertex we are looking at.
        int* arr = new int[count];    
  
        // earlier, I made every outofrange (technically 0 when
        // I first generated the matrix) value in my matrix
        // equal -1 so that it would be easy to work around.
        // This looks for every vertex attatched to the current value of i
        // and stores it in the dynamically allocated array.
        for (int j = 0; j < V - 1; j++)        
        {
            if (matrix[i][j] > -1)
            {
                //cout << matrix[i][j] << " ";
                arr[blob] = matrix[i][j];
                blob++;
            }           
        }
        // here I create the bounds of my array so that I can divide and find the mins and maxs.
        int R = count - 1;
        int L = 0;
        
        cout << "for " << e.A[i] << " " << endl;
        if (type == 'o')
        {
            // I need this condition to make sure I dont pass an empty array to my min/max functions
            if (count != 0)
            {
                cout << "the max out-degree is " << e.A[fun(arr, L, R, true)];
                cout << " and the min out-degree is " << e.A[fun(arr, L, R, false)] << endl;
            }
            else
                cout << e.A[i] << " does not have any out-degree vertices" << endl;
        }
        else 
        {
            if (count != 0)
            {
                cout << "the max in-degree is " << e.A[fun(arr, L, R, true)];
                cout << " and the min in-degree is " << e.A[fun(arr, L, R, false)] << endl;
            }
            else
                cout << e.A[i] << " does not have any in-degree vertices" << endl;
        
        }
        
    }
}

void printList(vector<int> adj[], int V, char type)
// Prints the list
{
    // I use an array of characters to make the output nice
    cout << "\n graph in list form" << endl;
    for (int d = 0; d < V - 1; ++d)
    {
        // depending on whether this is an undirected graph or directed graph,
        // the output will be different
        switch(type)
        {
        case 'u':
            // prints the active nodes            
            cout << "\n Vertex "
                << e.A[d] << ":";
            for (auto x : adj[d])
                cout << " " << e.A[x] << ",";
            printf("\n");
            break;
        case 'd':
            // prints the active nodes            
            cout << "\n Vertex "
                << e.A[d] << ":";
            for (auto x : adj[d])
                cout << " -> " << e.A[x];
            printf("\n");
            break;

        }
    }
}
int** Transposition(vector<int> adj[], int V)
// This takes a matrix of 0s and 1s and assigns
// its original values back into the 1s place.
// Everything else becomes -1.
{
    // declare the "2d array" and initialize its row
    int** matrix = new int* [V - 1];
    // iterates through to copy the 1s to my matrix
    for (int d = 0; d < V - 1; ++d)
    {
        // initialize its height
        matrix[d] = new int[V - 1];

        for (auto f : adj[d])
        {
            if (!adj[d].empty())
                matrix[d][f] = f;

        }
    }
    // Finishes up the matrix by assigning the rest as -1.
    // -1 is easier to work with than -853634563456345
    for (int i = 0; i < V - 1; i++)
    {
        for (int j = 0; j < V - 1; j++)
        {
            if (matrix[i][j] < -1)
                matrix[i][j] = -1;
        }
    }

    return matrix;
}
int** CreateMatrix(vector<int> adj[], int V)
// I soon had problems with an array of vectors since
// it wasn't the kindest thing to iterate through.
// That is why I copied over my data from the list (actually vector) to 
// a 2d array so I can have immediate access to my columns.
{
    // declare the "2d array" and initialize its row
    int** matrix = new int* [V - 1];
    // iterates through to copy the 1s to my matrix
    for (int d = 0; d < V - 1; ++d)
    {
        // initialize its height
        matrix[d] = new int[V - 1];

        for (auto f : adj[d])
        {
            if (!adj[d].empty())
                matrix[d][f] = 1;

        }
    }
    // Finishes up the matrix by assigning the rest as 0.
    // This solves my problem of iteration through my vectors
    // and alleviatses stupid outofbounds values.
    for (int i = 0; i < V - 1; i++)
    {
        for (int j = 0; j < V - 1; j++)
        {
            if (matrix[i][j] != 1)
                matrix[i][j] = 0;
        }
    }
    
    return matrix;
}
void printMatrix(int** matrix, int V) 
{
    // prints my matrix with row and column markers
        cout << "\n graph in matrix form" << endl;
        for (int i = 0; i < V - 1; i++)
            cout <<" "<< e.A[i] << " ";
        cout << "\n------------------  " << endl;
        for (int i = 0; i < V - 1; i++)
        {
            cout << endl;
            for (int j = 0; j < V - 1; j++)
            {
                cout << " " << matrix[i][j] << " ";
            }
            cout << " | " << e.A[i] << endl;
        }
        // the matrix should be symetric since it is undirected,
        // so if a node points at another node, the node its
        // pointing at (other node) will be pointing at the
        // node (first node).

        // If there isnt any self loops, then the diagnol should
        // be all zeros
}
void traverse(int** matrix, int currentVertex, bool visited[], int V) // V - 1 is the number of nodes that exist
// The next two functions are based on a Depth First Search.
{
    // We have found the current node!
    // This should prevent nodes already found to be 
    // discovered again since we are following
    // an arbitrary path
    visited[currentVertex] = true;

    // A queue that keeps track of what we found
    queue<int> found;

    // pushes the node to a queue
    found.push(currentVertex);

    // While we have nodes found, we continue 
    // to traverse the matrix 
    while (!found.empty())
    {
        // takes the next node and adds it the queue
        int nextVertex = found.front(); //delete from queue and print
        // kills the previous node
        found.pop();

        for (int i = 0; i < V - 1; i++)
        {
            // if the current value of the matrix is 1 (so true),
            // we ask if its been visited.
            if (matrix[i][nextVertex])
            {
                // if the node has not been visited, 
                // we mark it as visited and push it to the queue
                if (!visited[i])
                {
                    visited[i] = true;
                    found.push(i);
                }
            }
        }
    }
}

bool Connected(int** matrix, int V)
// this is the application of my DFS traversal function
{
    // we make a bool array to keep track of what we have visited
    bool* vis = new bool[V - 1];
    //starting at 0, it checks whether all nodes (vertices) are found or not
    for (int vertex = 0; vertex < V - 1; vertex++)
    {
        for (int i = 0; i < V - 1; i++)
        {
            // initialize as if a node has not been found. Primes it
            // for when we actually find a node.
            vis[i] = false;
        }
        
        traverse(matrix, vertex, vis, V);
        // pushes all the isolated nodes to the queue so that I can print them later
        for (int i = 0; i < V - 1; i++)
        {
             if (vis[i] == 0) 
                 e.isolatedNodes.push(i);
        }
        for (int i = 0; i < V - 1; i++)
        {
            // if there is a node not found by the traversal,
            // then the graph is not connected. Cuts off search when the first
            // isolated node is found;
            if (!vis[i]) 
            {                
                return false;
            }            
        }      
    }
    // otherwise, if all nodes are found the graph is connected
    return true;
}
void isConnected(int** matrix, int V)
{
    if (Connected(matrix, V))
        cout << "the graph is completely connected" << endl;
    else
        cout << "the graph is not completely connected" << endl;
}
void selfLoop(vector<int> adj[], int i)
// tells the user where self loops are
{
            cout << "There is a self loop at "
                 << e.A[i] << endl;
}
void isolatedNodes()
// tells the user where the isolated nodes are. 
// If none are there, then it outputs that there are none.
{
    if (!e.isolatedNodes.empty())
    {
        cout << "There is are isolated vertices at ";
        while (!e.isolatedNodes.empty()) {
            cout << e.A[e.isolatedNodes.front()] << ", ";
            e.isolatedNodes.pop();
        }
        cout << endl;
    }
    else
        cout << "There are no isolated nodes" << endl;

}
void checkEdges(vector<int> adj[], int V)
// checks the edges for self loops and tells the user if none exist
{
    // counts how many self loops are found
    int selfCount = 0;
    for (int i = 0; i < V - 1; i++)
        for (auto j : adj[i])
        {
            // increments if a self loop is found
            if (adj[i] == adj[j])
            {
                selfCount++;
                selfLoop(adj, i);       
            }
        }

    // checks if the graph has any self loops   
    if (selfCount == 0)
        cout << "\n The graph has no self loops" << endl;
    else
        cout << "\n There are a total of " << selfCount << " self loops" << endl;
}
void addEdge( int s, int d)
// Connects edges indirectly. Can be used directly if needed.
{
    switch (e.type) 
    {
    case 'u':
        // This can check for self loops, 
        // but the existance of self loops is 
        // not dependent on this. 
        if (s != d)
        {
            // connects node s to node d
            e.adj[s].push_back(d);
            // connects node d to node s
            e.adj[d].push_back(s);
            // thus connects in both 
            // directions as an undirected graph       
        }
        else
        {
            // connects node s to node d
            e.adj[s].push_back(d);
            // thus connects in one
            // direction as a directed graph 

            // This logic is also used for self loops
            // so that duplicates do not poison the vector 
        }  
        break;
    case 'd':
        e.adj[s].push_back(d);
        e.trans[d].push_back(s);
        break;
    }
    // increment the edge counter by 1.
    e.edgeNum++;
}
int main()
{
    cout << "is your graph direct 'd' or undirect 'u'" << endl;
    cin >> e.type;

    // initializes the edge counter to 0
    e.edgeNum = 0;
    e.doubleCount = 0;

    // Size of vector. Number of possible nodes = V - 1. 
    // For simplicities sake, we will limit the 
    // vertices to the uppercase letters of the alphabet.
    int V = 7;

    // adding all of the edges to the graph    
    addEdge(0, 1);
    addEdge(1, 0); // This is used to show that A and B are mutually adjacent
    //addEdge(0, 0); // A has a self loop
    addEdge(0, 3);
    addEdge(2, 2);
    addEdge(3, 2);
    //addEdge(3, 1);
    //addEdge(1, 4);
    addEdge(4, 1);
    addEdge(0, 5);  // F is now a isolated vertex

    // declaration and initialization of my matrices
    //----------------------------------------------
    // this matrix contains 1s and 0s representing the adjency list
    int** adjMatrix = CreateMatrix(e.adj, V);
    // this is a transpose of adjmatrix 
    int** transMatrix = CreateMatrix(e.trans, V);

    // This matrix contains its original values of vector adj
    // but replaces outofbounds values with -1s.
    int** reverseAdj= Transposition(e.adj,  V);
    // Similar to reversematrix, it does the same thing, but for the reversed adj 
    int** reverseTrans = Transposition(e.trans, V);

    //results------------------------------------------
    switch (e.type)
    {
    case 'u':
        // prints the adjency list
        printList(e.adj, V, e.type);
        // checks the edges for self loops and tells the user if none exist
        checkEdges(e.adj, V);
        // checks if the graph is connected. 
        isConnected(adjMatrix, V);
        // displays what nodes are isolated
        isolatedNodes();
        break;
    case 'd':
        // prints the adjency matrix
        printMatrix(adjMatrix, V);
        // displays out-degree for each vertex
        Degree(reverseAdj, V, 'o');
        // displays in-deg for each vertex
        Degree(reverseTrans, V, 'i');
        adjacentVertex(adjMatrix, transMatrix, V);
        break;
    }

    cout << "number of edges " << e.edgeNum << endl;

    cout << "There are a total of " << V - 1 << " active nodes" << endl;


}
