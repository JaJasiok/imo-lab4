#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <climits>
#include <algorithm>
#include <unordered_set>
#include <random>
#include <numeric>
#include <variant>
#include <chrono>
#include <unistd.h>

using namespace std;

vector<vector<int>> readKroaFile(const string &filename)
{
    ifstream file(filename);
    string line;
    vector<vector<int>> verticesCoords;

    if (file.is_open())
    {
        // Skip the header lines
        for (int i = 0; i < 7; i++)
        {
            getline(file, line);
        }

        // Read the coordinates and populate the cost matrix
        do
        {
            istringstream iss(line);
            int index, x, y;
            iss >> index >> x >> y;
            verticesCoords.push_back({x, y});
            // cout << index << " " << x << " " << y << endl;
        } while (getline(file, line) && line != "EOF");

        // verticesCoords.pop_back();

        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }

    return verticesCoords;
}

vector<vector<int>> createDistanceMatrix(const vector<vector<int>> &verticesCoords)
{
    int numVertices = verticesCoords.size();
    vector<vector<int>> distanceMatrix(numVertices, vector<int>(numVertices));

    for (int i = 0; i < numVertices; i++)
    {
        for (int j = 0; j < numVertices; j++)
        {
            double x1 = verticesCoords[i][0];
            double y1 = verticesCoords[i][1];
            double x2 = verticesCoords[j][0];
            double y2 = verticesCoords[j][1];

            double distance = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
            distanceMatrix[i][j] = int(round(distance));

            // cout << i << " " << j << " " << distanceMatrix[i][j] << endl;
        }
    }

    return distanceMatrix;
}

class Vertex
{
public:
    int id;
    Vertex *prev;
    Vertex *next;

    Vertex(int _id) : id(_id), prev(nullptr), next(nullptr) {}
};

class Edge
{
public:
    Vertex *src;
    Vertex *dest;
    int distance;

    Edge(Vertex *_src, Vertex *_dest, int _distance) : src(_src), dest(_dest), distance(_distance) {}

    void remove()
    {
        delete this;
    }
};

class Graph
{
public:
    vector<Vertex *> vertices;
    vector<Edge *> edges;
    int distance = 0;

    Graph() {}

    Graph(const Graph& other)
    {
        for (Vertex* v : other.vertices)
        {
            Vertex* newVertex = new Vertex(v->id);
            addVertex(newVertex);
        }

        for (Edge* e : other.edges)
        {
            Vertex* src = findVertex(e->src->id);
            Vertex* dest = findVertex(e->dest->id);
            addEdge(src, dest, e->distance);
        }
    }

    void addVertex(Vertex *v)
    {
        if (v)
        {
            if (find(vertices.begin(), vertices.end(), v) == vertices.end())
            {
                vertices.push_back(v);
            }
        }
    }

    void addEdge(Vertex *src, Vertex *dest, int distance = 0)
    {
        try
        {
            if (src && dest && src != dest)
            {
                for (Edge *e : edges)
                {
                    if ((e->src == src && e->dest == dest) || (e->src == dest && e->dest == src))
                    {
                        return;
                    }
                }

                Edge *e = new Edge(src, dest, distance);
                src->next = dest;
                dest->prev = src;
                this->distance += distance;
                edges.push_back(e);
            }
        }
        catch (exception &e)
        {
            cout << e.what() << endl;
        }
    }

    void removeVertex(Vertex *v)
    {
        if (v)
        {
            removeEdge(v->prev, v);
            removeEdge(v, v->next);

            vertices.erase(remove(vertices.begin(), vertices.end(), v), vertices.end());
        }
    }

    void removeEdge(Vertex *src, Vertex *dest)
    {
        if (src && dest)
        {
            for (Edge *e : edges)
            {
                if (e->src == src && e->dest == dest)
                {
                    this->distance -= e->distance;
                    edges.erase(remove(edges.begin(), edges.end(), e), edges.end());
                    delete e;
                    break;
                }
            }
            src->next = nullptr;
            dest->prev = nullptr;
        }
    }

    Vertex *findVertex(int id)
    {
        for (Vertex *v : vertices)
        {
            if (v->id == id)
                return v;
        }
        return nullptr;
    }

    Edge *findEdge(int src_id, int dest_id)
    {
        for (Edge *e : edges)
        {
            if (e->src->id == src_id && e->dest->id == dest_id)
            {
                return e;
            }
        }
        return nullptr;
    }

    vector<Edge *> findPath(Vertex *start, Vertex *end)
    {
        vector<Edge *> path;
        Vertex *currentVertex = start;

        do
        {
            for (Edge *e : edges)
            {
                if (e->src == currentVertex)
                {
                    path.push_back(e);
                    currentVertex = e->dest;
                    break;
                }
            }
        } while (currentVertex != end);

        return path;
    }

    void destroy(Vertex *start, int n, const vector<vector<int>> &distanceMatrix){
        Vertex *currentVertex = start;
        for(int i = 0; i < n/2; i++){
            currentVertex = currentVertex->prev;
        }
        Vertex *startPrev = currentVertex->prev;
        for(int i = 0; i < n; i++){
            Vertex *next = currentVertex->next;

            removeVertex(currentVertex);
            currentVertex = next;
            if(i == n - 1){
                addEdge(startPrev, currentVertex, distanceMatrix[startPrev->id][currentVertex->id]);
            }
            // if(i != 0){
            //     delete currentVertex->prev;
            // }
        }
    }
};

class Move
{
public:
    pair<Vertex *, Vertex *> vertices;
    pair<Edge *, Edge *> edges;
    Graph *graph;
    int delta;

    Move(pair<Vertex *, Vertex *> _vertices = {nullptr, nullptr}, pair<Edge *, Edge *> _edges = {nullptr, nullptr}, Graph *_graph = nullptr, int _delta = 0)
        : vertices(_vertices), edges(_edges), graph(_graph), delta(_delta) {}
};

void saveGraphs(const vector<Graph> &graphs, const string &filename)
{
    ofstream file(filename);
    if (file.is_open())
    {
        for (Graph g : graphs)
        {
            for (Edge *e : g.edges)
            {
                file << e->src->id << " " << e->dest->id << endl;
            }
            file << endl;
        }
        file.close();
    }
    else
    {
        cout << "Failed to open file: " << filename << endl;
    }
}

vector<Graph> randomCycles(const vector<vector<int>> &distanceMatrix)
{
    int numVertices = distanceMatrix.size();

    std::vector<int> numbers(numVertices);
    iota(numbers.begin(), numbers.end(), 0);

    random_shuffle(numbers.begin(), numbers.end());

    vector<Graph> cycles(2);

    std::vector<int> group1(numbers.begin(), numbers.begin() + numVertices / 2);
    std::vector<int> group2(numbers.begin() + numVertices / 2, numbers.end());

    for (int i : group1)
    {
        Vertex *v = new Vertex(i);
        cycles[0].addVertex(v);
    }

    for (int i : group2)
    {
        Vertex *v = new Vertex(i);
        cycles[1].addVertex(v);
    }

    for (int i = 0; i < numVertices / 2; i++)
    {
        if (i < (numVertices / 2) - 1)
        {
            Vertex *v1 = cycles[0].vertices[i];
            Vertex *v2 = cycles[0].vertices[i + 1];
            // cout << v1->id << " " << v2->id << " " << distanceMatrix[v1->id][v2->id] << endl;
            cycles[0].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
        else
        {
            Vertex *v1 = cycles[0].vertices[i];
            Vertex *v2 = cycles[0].vertices[0];
            // cout << v1->id << " " << v2->id << " " << distanceMatrix[v1->id][v2->id] << endl;
            cycles[0].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
    }

    for (int i = 0; i < numVertices / 2; i++)
    {
        if (i < (numVertices / 2) - 1)
        {
            Vertex *v1 = cycles[1].vertices[i];
            Vertex *v2 = cycles[1].vertices[i + 1];
            // cout << distanceMatrix[v1->id][v2->id] << endl;
            cycles[1].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
        else
        {
            Vertex *v1 = cycles[1].vertices[i];
            Vertex *v2 = cycles[1].vertices[0];
            // cout << distanceMatrix[v1->id][v2->id] << endl;
            cycles[1].addEdge(v1, v2, distanceMatrix[v1->id][v2->id]);
        }
    }

    return cycles;
}

void swapVerticesBetweenCycles(Vertex *vertex1, Vertex *vertex2, vector<Graph> &graphs, const vector<vector<int>> &distanceMatrix)
{
    Graph *graph1;
    Graph *graph2;

    if(graphs[0].findVertex(vertex1->id) != nullptr)
    {
        graph1 = &graphs[0];
        graph2 = &graphs[1];
    }
    else
    {
        graph1 = &graphs[1];
        graph2 = &graphs[0];
    }


    Vertex *vertex1Prev = vertex1->prev;
    Vertex *vertex1Next = vertex1->next;

    Vertex *vertex2Prev = vertex2->prev;
    Vertex *vertex2Next = vertex2->next;

    graph1->removeVertex(vertex1);

    graph2->removeVertex(vertex2);

    graph1->addVertex(vertex2);
    graph1->addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
    graph1->addEdge(vertex2, vertex1Next, distanceMatrix[vertex2->id][vertex1Next->id]);

    graph2->addVertex(vertex1);
    graph2->addEdge(vertex2Prev, vertex1, distanceMatrix[vertex2Prev->id][vertex1->id]);
    graph2->addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
}

void swapEdgesInGraph(Edge *edge1, Edge *edge2, Graph *graph, const vector<vector<int>> &distanceMatrix)
{
    Vertex *vertex11 = edge1->src;
    Vertex *vertex12 = edge1->dest;
    Vertex *vertex21 = edge2->src;
    Vertex *vertex22 = edge2->dest;

    // for(Edge *e: graph->edges)
    // {
    //     cout << e->src->id << " " << e->dest->id << endl;
    // }

    graph->removeEdge(vertex11, vertex12);
    graph->removeEdge(vertex21, vertex22);

    vector<Edge *> pathToReverse = graph->findPath(vertex12, vertex21);

    for (Edge *edge : pathToReverse)
    {
        swap(edge->src, edge->dest);
        swap(edge->dest->next, edge->dest->prev);
        if (edge->src->next == nullptr)
        {
            swap(edge->src->next, edge->src->prev);
        }
    }

    graph->addEdge(vertex11, vertex21, distanceMatrix[vertex11->id][vertex21->id]);
    graph->addEdge(vertex12, vertex22, distanceMatrix[vertex12->id][vertex22->id]);
}

void swapVerticesInCycle(Vertex *vertex1, Vertex *vertex2, Graph *graph, const vector<vector<int>> &distanceMatrix)
{
    Vertex *vertex1Prev = vertex1->prev;
    Vertex *vertex1Next = vertex1->next;

    Vertex *vertex2Prev = vertex2->prev;
    Vertex *vertex2Next = vertex2->next;


    graph->removeEdge(vertex1Prev, vertex1);
    graph->removeEdge(vertex1, vertex1Next);
    graph->removeEdge(vertex2Prev, vertex2);
    graph->removeEdge(vertex2, vertex2Next);

    graph->addEdge(vertex1Prev, vertex2, distanceMatrix[vertex1Prev->id][vertex2->id]);
    graph->addEdge(vertex2, vertex1Next, distanceMatrix[vertex2->id][vertex1Next->id]);
    graph->addEdge(vertex2Prev, vertex1, distanceMatrix[vertex2Prev->id][vertex1->id]);
    graph->addEdge(vertex1, vertex2Next, distanceMatrix[vertex1->id][vertex2Next->id]);
}

void repair(vector<Graph>& cycles, const vector<vector<int>> &distanceMatrix)
{
    int numVertices = distanceMatrix.size();
    vector<bool> visited(numVertices, false);

    for (Graph &cycle : cycles)
    {
        for (Vertex *v : cycle.vertices)
        {
            visited[v->id] = true;
        }
    }

    do
    {
        int minDistance = INT_MAX;
        int vertexId = -1;
        pair<Edge *, Graph *> minPair;

        vector<pair<Edge *, Graph *>> edgesInGraphs;
        if (cycles[0].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[1].edges)
            {
                edgesInGraphs.push_back({e, &cycles[1]});
            }
        }
        else if (cycles[1].vertices.size() == numVertices / 2)
        {
            for (Edge *e : cycles[0].edges)
            {
                edgesInGraphs.push_back({e, &cycles[0]});
            }
        }
        else
        {
            for (Graph &cycle : cycles)
            {
                for (Edge *e : cycle.edges)
                {
                    edgesInGraphs.push_back({e, &cycle});
                }
            }
        }

        for (int j = 0; j < numVertices; j++)
        {
            if (!visited[j])
            {
                for (pair<Edge *, Graph *> singlePair : edgesInGraphs)
                {
                    Edge *e = singlePair.first;
                    int distanceSum = cycles[0].distance + cycles[1].distance;
                    if (distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j] < minDistance)
                    {
                        minDistance = distanceSum - e->distance + distanceMatrix[e->dest->id][j] + distanceMatrix[e->src->id][j];
                        vertexId = j;
                        minPair = singlePair;
                    }
                }
            }
        }

        Graph *cyclePtr = minPair.second;
        Graph &cycle = *cyclePtr;

        Edge *minEdge = minPair.first;

        Vertex *newVertex = new Vertex(vertexId);

        cycle.addVertex(newVertex);
        cycle.addEdge(minEdge->src, newVertex, distanceMatrix[minEdge->src->id][vertexId]);
        cycle.addEdge(newVertex, minEdge->dest, distanceMatrix[newVertex->id][minEdge->dest->id]);

        if (cycle.vertices.size() > 3)
        {
            cycle.removeEdge(minEdge->src, minEdge->dest);
        }

        visited[vertexId] = true;
    }while(cycles[0].vertices.size() < numVertices / 2 || cycles[1].vertices.size() < numVertices / 2);
}

vector<Graph> steepestLocalSearch(const vector<Graph> &startCycles, const vector<vector<int>> &distanceMatrix)
{
    vector<Graph> cycles;

    for(Graph g: startCycles)
    {
        Graph newGraph(g);
        cycles.push_back(Graph(g));
    }


    int bestDelta = 0;

    int i = 0;

    do
    {
        vector<Move *> moves;

        for (Vertex *vertex1 : cycles[0].vertices)
        {
            for (Vertex *vertex2 : cycles[1].vertices)
            {
                int delta = distanceMatrix[vertex1->prev->id][vertex2->id] + distanceMatrix[vertex2->id][vertex1->next->id] + distanceMatrix[vertex2->prev->id][vertex1->id] + distanceMatrix[vertex1->id][vertex2->next->id] - distanceMatrix[vertex1->prev->id][vertex1->id] - distanceMatrix[vertex1->id][vertex1->next->id] - distanceMatrix[vertex2->prev->id][vertex2->id] - distanceMatrix[vertex2->id][vertex2->next->id];

                if (delta < 0)
                {
                    Move *move = new Move({vertex1, vertex2}, {nullptr, nullptr}, nullptr, delta);
                    moves.push_back(move);
                }
            }
        }

        for (Graph &graph : cycles)
        {
            for (Edge *edge1 : graph.edges)
            {
                for (Edge *edge2 : graph.edges)
                {
                    if (edge1 != edge2)
                    {
                        Vertex *vertex11 = edge1->src;
                        Vertex *vertex12 = edge1->dest;
                        Vertex *vertex21 = edge2->src;
                        Vertex *vertex22 = edge2->dest;

                        if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                        {
                            int delta = distanceMatrix[vertex11->id][vertex21->id] + distanceMatrix[vertex12->id][vertex22->id] - edge1->distance - edge2->distance;
                            if (delta < 0)
                            {
                                Move *move = new Move({nullptr, nullptr}, {edge1, edge2}, &graph, delta);
                                moves.push_back(move);
                            }
                        }
                    }
                }
            }
        }

        if(moves.size() > 0){

            sort(moves.begin(), moves.end(), [](Move *a, Move *b)
                { return a->delta < b->delta; });

            Move *bestMove = moves[0];

            if (bestMove->graph == nullptr)
            {
                Vertex *vertex1 = bestMove->vertices.first;
                Vertex *vertex2 = bestMove->vertices.second;

                swapVerticesBetweenCycles(vertex1, vertex2, cycles, distanceMatrix);
            }
            else
            {
                Edge *edge1 = bestMove->edges.first;
                Edge *edge2 = bestMove->edges.second;

                Graph *graphToSwap = bestMove->graph;

                swapEdgesInGraph(edge1, edge2, graphToSwap, distanceMatrix);
            }

            for (Move* move : moves) {
                delete move;
            }
            moves.clear();
        }
        else{
            break;
        }

        i++;
    } while (true);

    return cycles;
}

vector<Graph> multipleLocalSearch(const vector<vector<int>> &distanceMatrix)
{
    vector<Graph> bestCycles;
    int bestDistance = INT_MAX;

    for (int i = 0; i < 100; i++)
    {
        vector<Graph> cycles = randomCycles(distanceMatrix);

        cycles = steepestLocalSearch(cycles, distanceMatrix);

        int newDistance = cycles[0].distance + cycles[1].distance;

        if (newDistance < bestDistance)
        {
            bestDistance = newDistance;
            bestCycles.clear();
            for(Graph g: cycles)
            {
                Graph newGraph(g);
                bestCycles.push_back(Graph(g));
            }
        }
    }

    return bestCycles;
}

vector<Graph> iteratedLocalSearch1(const vector<Graph> &startCycles, const vector<vector<int>> &distanceMatrix, int timeLimit)
{
    vector<Graph> bestCycles;

    for(Graph g: startCycles)
    {
        Graph newGraph(g);
        bestCycles.push_back(Graph(g));
    }

    bestCycles = steepestLocalSearch(bestCycles, distanceMatrix);

    int a = 0;

    int elapsed;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    do{
        vector<Graph> newCycles;

        for(Graph g: bestCycles)
        {
            Graph newGraph(g);
            newCycles.push_back(Graph(g));
        }

        // Distance: 34139 - 10
        // Distance: 32975 - 15
        // Distance: 33079 - 20

        for(int i = 0; i < 15; i ++){

            vector<Move *> moves;

            for (Vertex *vertex1 : newCycles[0].vertices)
            {
                for (Vertex *vertex2 : newCycles[1].vertices)
                {
                    Move *move = new Move({vertex1, vertex2}, {nullptr, nullptr}, nullptr, 0);
                    moves.push_back(move);
                }
            }

            for (Graph &graph : newCycles)
            {
                for (Edge *edge1 : graph.edges)
                {
                    for (Edge *edge2 : graph.edges)
                    {
                        if (edge1 != edge2)
                        {
                            Vertex *vertex11 = edge1->src;
                            Vertex *vertex12 = edge1->dest;
                            Vertex *vertex21 = edge2->src;
                            Vertex *vertex22 = edge2->dest;

                            if (vertex11->id != vertex12->id && vertex11->id != vertex21->id && vertex11->id != vertex22->id && vertex12->id != vertex21->id && vertex12->id != vertex22->id && vertex21->id != vertex22->id)
                            {
                                Move *move = new Move({nullptr, nullptr}, {edge1, edge2}, &graph, 0);
                                moves.push_back(move);
                            }
                        }
                    }
                }
            }

            Move* randomMove = moves[rand() % moves.size()];

            if (randomMove->graph == nullptr)
            {
                Vertex *vertex1 = randomMove->vertices.first;
                Vertex *vertex2 = randomMove->vertices.second;

                swapVerticesBetweenCycles(vertex1, vertex2, newCycles, distanceMatrix);
            }
            else
            {
                Edge *edge1 = randomMove->edges.first;
                Edge *edge2 = randomMove->edges.second;

                Graph *graphToSwap = randomMove->graph;

                swapEdgesInGraph(edge1, edge2, graphToSwap, distanceMatrix);
            }

            for (Move* move : moves) {
                delete move;
            }
            moves.clear();
        }

        newCycles = steepestLocalSearch(newCycles, distanceMatrix);

        if(newCycles[0].distance + newCycles[1].distance < bestCycles[0].distance + bestCycles[1].distance)
        {
            bestCycles.clear();
            for(Graph g: newCycles)
            {
                Graph newGraph(g);
                bestCycles.push_back(Graph(g));
            }
        }

        a++;

        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
    } while(elapsed < timeLimit);

    cout << a << endl;

    return bestCycles;
}

vector<Graph> iteratedLocalSearch2a(const vector<Graph> &startCycles, const vector<vector<int>> &distanceMatrix, int timeLimit)
{
    vector<Graph> bestCycles;

    for(Graph g: startCycles)
    {
        Graph newGraph(g);
        bestCycles.push_back(Graph(g));
    }

    int a = 0;

    int elapsed;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    bestCycles = steepestLocalSearch(bestCycles, distanceMatrix);

    do{
        vector<Graph> newCycles;

        for(Graph g: bestCycles)
        {
            Graph newGraph(g);
            newCycles.push_back(Graph(g));
        }

        // int minDistance = INT_MAX;

        // Vertex *nearestVertex1 = nullptr;
        // Vertex *nearestVertex2 = nullptr;

        // for(Vertex *vertex1: newCycles[0].vertices)
        // {
        //     for(Vertex *vertex2: newCycles[1].vertices)
        //     {
        //         if(distanceMatrix[vertex1->id][vertex2->id] < minDistance)
        //         {
        //             minDistance = distanceMatrix[vertex1->id][vertex2->id];
        //             nearestVertex1 = vertex1;
        //             nearestVertex2 = vertex2;
        //         }
        //     }
        // }

        newCycles[0].destroy(newCycles[0].vertices[rand() % newCycles[0].vertices.size()], 30, distanceMatrix);
        newCycles[1].destroy(newCycles[1].vertices[rand() % newCycles[1].vertices.size()], 30, distanceMatrix);

        repair(newCycles, distanceMatrix);

        if(newCycles[0].distance + newCycles[1].distance < bestCycles[0].distance + bestCycles[1].distance)
        {
            bestCycles.clear();
            for(Graph g: newCycles)
            {
                Graph newGraph(g);
                bestCycles.push_back(Graph(g));
            }
        }

        a++;

        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
    } while(elapsed < timeLimit);

    cout << a << endl;

    return bestCycles;
}

vector<Graph> iteratedLocalSearch2b(const vector<Graph> &startCycles, const vector<vector<int>> &distanceMatrix, int timeLimit)
{
    vector<Graph> bestCycles;

    for(Graph g: startCycles)
    {
        Graph newGraph(g);
        bestCycles.push_back(Graph(g));
    }

    int a = 0;

    int elapsed;

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    do{
        vector<Graph> newCycles;

        for(Graph g: bestCycles)
        {
            Graph newGraph(g);
            newCycles.push_back(Graph(g));
        }

        newCycles = steepestLocalSearch(newCycles, distanceMatrix);

        // int minDistance = INT_MAX;

        // Vertex *nearestVertex1 = nullptr;
        // Vertex *nearestVertex2 = nullptr;

        // for(Vertex *vertex1: newCycles[0].vertices)
        // {
        //     for(Vertex *vertex2: newCycles[1].vertices)
        //     {
        //         if(distanceMatrix[vertex1->id][vertex2->id] < minDistance)
        //         {
        //             minDistance = distanceMatrix[vertex1->id][vertex2->id];
        //             nearestVertex1 = vertex1;
        //             nearestVertex2 = vertex2;
        //         }
        //     }
        // }

        newCycles[0].destroy(newCycles[0].vertices[rand() % newCycles[0].vertices.size()], 30, distanceMatrix);
        newCycles[1].destroy(newCycles[1].vertices[rand() % newCycles[1].vertices.size()], 30, distanceMatrix);

        repair(newCycles, distanceMatrix);

        if(newCycles[0].distance + newCycles[1].distance < bestCycles[0].distance + bestCycles[1].distance)
        {
            bestCycles.clear();
            for(Graph g: newCycles)
            {
                Graph newGraph(g);
                bestCycles.push_back(Graph(g));
            }
        }

        a++;

        chrono::steady_clock::time_point end = chrono::steady_clock::now();

        elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();
    } while(elapsed < timeLimit);

    cout << a << endl;

    return bestCycles;
}

int main()
{
    // srand(time(0));

    vector<vector<int>> verticesCoords = readKroaFile("kroA200.tsp");

    vector<vector<int>> distanceMatrix = createDistanceMatrix(verticesCoords);

    vector<Graph> startCycles = randomCycles(distanceMatrix);

    chrono::steady_clock::time_point begin = chrono::steady_clock::now();

    vector<Graph> multi = multipleLocalSearch(distanceMatrix);

    chrono::steady_clock::time_point end = chrono::steady_clock::now();

    int elapsed = chrono::duration_cast<chrono::milliseconds>(end - begin).count();

    cout << elapsed << endl;

    cout << "Distance: " << multi[0].distance + multi[1].distance << endl;

    vector<Graph> iterated1 = iteratedLocalSearch1(startCycles, distanceMatrix, elapsed);

    cout << "Distance: " << iterated1[0].distance + iterated1[1].distance << endl;

    vector<Graph> iterated2a = iteratedLocalSearch2a(startCycles, distanceMatrix, elapsed);

    cout << "Distance: " << iterated2a[0].distance + iterated2a[1].distance << endl;

    vector<Graph> iterated2b = iteratedLocalSearch2b(startCycles, distanceMatrix, elapsed);

    cout << "Distance: " << iterated2b[0].distance + iterated2b[1].distance << endl;

    return 0;
}