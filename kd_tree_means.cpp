#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>  
#include <random> 
#include <iomanip> 
using namespace std;

vector<double*> points;

vector<double*> load_points_from_csv(string filename, int dimension) {
    ifstream file(filename);
    vector<double*> points;

    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo " << filename << endl;
        return points;
    }

    string line;

    while (getline(file, line)) {
        istringstream ss(line);
        double* point = new double[dimension];

        for (int i = 0; i < dimension; ++i) {
            string value;
            if (!getline(ss, value, ',')) break;
            point[i] = stod(value);
        }
        points.push_back(point);
    }

    file.close();
    return points;
}

class KDtree {


public:
    int dimension;
    struct Node {
        Node* left;
        Node* right;
        double* point;
        Node(int dimension) : left(nullptr), right(nullptr) { point = new double[dimension]; }
    };
    Node* root;

    KDtree(int dimension) : dimension(dimension), root(nullptr) {}

    void insert(double* point) {
        root = insert_node(root, point, 0);
    }

    double* nearest(double* point) {
        double best_dist = 10000000;
        return search_near(root, point, 0, best_dist);
    }

    double euclidian_dist(double* p1, double* p2) {
        double sum = 0;
        for (int i = 0; i < dimension; i++) {
            sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
        }
        return sqrt(sum);
    }

    void print_node(Node* node, int depth) {
        if (node == nullptr) return;

        std::cout << std::setw(depth * 2) << ""; // Indentación según profundidad
        std::cout << "Nivel " << depth << " | Punto: (";
        for (int i = 0; i < dimension; i++) {
            std::cout << node->point[i];
            if (i < dimension - 1) std::cout << ", ";
        }
        std::cout << ")" << std::endl;

        // Recursivamente imprimir los hijos
        if (node->left != nullptr) {
            std::cout << std::setw(depth * 2) << "" << "Izquierda -> ";
            print_node(node->left, depth + 1);
        }
        if (node->right != nullptr) {
            std::cout << std::setw(depth * 2) << "" << "Derecha -> ";
            print_node(node->right, depth + 1);
        }
    }

private:

    Node* insert_node(Node* node, double* point, int dim) {
        if (node == nullptr) {
            Node* temp = new Node(dimension);
            for (int i = 0; i < dimension; i++) {
                temp->point[i] = point[i];
            }
            return temp;
        }
        int div = dim % dimension;
        if (point[div] < node->point[div]) {
            node->left = insert_node(node->left, point, dim + 1);
        }
        else {
            node->right = insert_node(node->right, point, dim + 1);
        }
        return node;
    }

    double* search_near(Node* node, double* point, int dim, double& best_dist) {
        if (node == nullptr) {
            return nullptr;
        }
        double dist = euclidian_dist(node->point, point);
      //  cout << node->point[0] << " - " << node->point[1] << "||" << point[0] << " - " << point[1] << "||" << dist << endl;
        double* best_point = nullptr;
        if (dist < best_dist) {
            best_dist = dist;
            best_point = node->point;
        }
        int div = dim % dimension;
        Node* prolly_branch = (point[div] < node->point[div]) ? node->left : node->right;
        Node* other_branch = (point[div] < node->point[div]) ? node->right : node->left;

        double* next_best_point = search_near(prolly_branch, point, dim + 1, best_dist);
        if (next_best_point != nullptr) {
            best_point = next_best_point;
        }
        if (abs(point[div] - node->point[div]) < best_dist) {
            double* other_best_point = search_near(other_branch, point, dim + 1, best_dist);
            if (other_best_point != nullptr && euclidian_dist(other_best_point, point) < best_dist) {
                best_dist = euclidian_dist(other_best_point, point);
                best_point = other_best_point;
            }
        }

        return best_point;
    }





};

class Kmeans {
    int k;
    int dimension;
    KDtree centroid_tree;

public:
    Kmeans(int k, int dimension)
        : k(k), dimension(dimension), centroid_tree(dimension) {
        srand(time(0));
    }

    vector<double*> centroids;

    void initialize() {
        centroids.clear();
        for (int i = 0; i < k; ++i) {
            centroids.push_back(points[rand() % points.size()]);
            centroid_tree.insert(centroids.back());
        }
        centroid_tree.print_node(centroid_tree.root, 0);
    }

    double* closest_centroid(double* point) {
        return centroid_tree.nearest(point);
    }

    void update_centroids(const vector<vector<double*>>& clusters) {
        vector<double*> new_centroids;
        for (int i = 0; i < k; ++i) {
            double* new_centroid = new double[dimension]();
            int cluster_size = clusters[i].size();

            if (cluster_size == 0) {
                cerr << "Advertencia: Cluster " << i << " está vacío. Asignando un punto aleatorio." << endl;
                new_centroid = points[rand() % points.size()];
            }
            else {
                for (auto& point : clusters[i]) {
                    for (int d = 0; d < dimension; ++d) {
                        new_centroid[d] += point[d];
                    }
                }
                for (int d = 0; d < dimension; ++d) {
                    new_centroid[d] /= cluster_size;
                }
            }
            new_centroids.push_back(new_centroid);
        }

        centroids = new_centroids;
        centroid_tree = KDtree(dimension);
        for (auto& centroid : centroids) {
            centroid_tree.insert(centroid);
        }
    }

    void fit() {
        initialize();
        bool converged = false;
        vector<vector<double*>> clusters(k);

        while (!converged) {
            for (int i = 0; i < k; ++i) {
                clusters[i].clear();
            }
            for (auto& point : points) {
                double* nearest_centroid = closest_centroid(point);
                for (int i = 0; i < k; ++i) {
                    if (euclidian_equal(nearest_centroid, centroids[i])) {
                        clusters[i].push_back(point);
                        break;
                    }
                }
            }

            vector<double*> old_centroids = centroids;
            update_centroids(clusters);
            converged = true;

            for (int i = 0; i < k; ++i) {
                if (centroid_tree.euclidian_dist(old_centroids[i], centroids[i]) > 0.0001) {
                    converged = false;
                    break;
                }
            }
        }

        save_clusters_to_csv(clusters);
    }

    void save_clusters_to_csv(const vector<vector<double*>>& clusters) {
        ofstream file("clusters.csv");
        file << "x,y,cluster_id" << endl;

        for (int cluster_id = 0; cluster_id < clusters.size(); ++cluster_id) {
            for (auto& point : clusters[cluster_id]) {
                file << point[0] << "," << point[1] << "," << cluster_id << endl;
            }
        }

        file.close();
    }

    void print_centroids() {
        centroid_tree.print_node(centroid_tree.root, 0);
    }

private:
    bool euclidian_equal(double* p1, double* p2) {
        return centroid_tree.euclidian_dist(p1, p2) < 1e-6;
    }
};



int main() {
    int dimension = 2;

    int k = 150; // Num de centroides 
    points = load_points_from_csv("data2k.csv", dimension);
    Kmeans kmeans(k, dimension);
    kmeans.fit();
    //kmeans.print_centroids();

    return 0;
}