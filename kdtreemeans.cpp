#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <ctime>
#include <map>
#include <iomanip> 
using namespace std;

vector<double*> load_points_from_csv(string& filename, int dimension) {
    ifstream file(filename);
    vector<double*> points;
    
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo " << filename << endl;
        return points;
    }
    
    string line;
    getline(file, line); // Ignorar encabezado
    
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

class KDtree{
   

    public:
    int dimension;
    struct Node{
        Node* left;
        Node* right;
        double* point;
        Node(int dimension): left(nullptr), right(nullptr){point = new double[dimension];}
    };
    Node* root;

    KDtree(int dimension) : dimension(dimension), root(nullptr) {}

    void insert(double* point){
        root=insert_node(root,point,0);
    }

    double* nearest(double* point) {
        double best_dist = 10000000;
        return search_near(root, point, 0, best_dist);
    }
        
    double euclidian_dist(double* p1, double* p2){
        double sum=0;
        for(int i=0;i<dimension;i++){
            sum+=(p1[i]-p2[i])*(p1[i]-p2[i]);
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

    Node* insert_node(Node* node,double* point,int dim){
        if(node==nullptr){
            Node* temp= new Node (dimension);
            for(int i=0;i<dimension;i++){
                temp->point[i]=point[i];
            }
            return temp;
        }

        int div=dim%dimension; // se refiere a en que dimension se hara la division teniendo en cuenta el total de dimensiones del tree

        if(point[div]<node->point[div]){
            node->left=insert_node(node->left,point,dim+1);
        }
        else{
            node->right=insert_node(node->right,point,dim+1);
        }
        return node;
    }

    double* search_near(Node* node,double* point, int dim,double& best_dist){
        if(node==nullptr){return nullptr;}
        double dist=euclidian_dist(node->point,point);
        double* best_point=nullptr;

        if(dist<best_dist){best_dist=dist; best_point=node->point;}

        int div=dim%dimension;
        Node* prolly_branch=(point[div]<node->point[div]) ? node->left:node->right;
        Node* other_branch=(point[div]<node->point[div]) ? node->right:node->left;

        double* next_best_point=search_near(prolly_branch,point,dim+1,best_dist);

        if(next_best_point!=nullptr&&euclidian_dist(next_best_point, point) < best_dist){best_point=next_best_point;}
        
        if(abs(point[div]-node->point[div])<best_dist){
            double* other_best_point=search_near(other_branch,point,dim+1,best_dist);
            if (other_best_point != nullptr&&euclidian_dist(other_best_point, point) < best_dist) {
                best_point = other_best_point;
            }
        }
        return best_point;
    }

    


};

class Kmeans{
    int k; 
    int dimension; 
    KDtree centroid_tree; 
    vector<double*> centroids; 
    map<double*, int> centroid_to_cluster_id; 
    public:
    Kmeans(int k, int dimension):k(k), dimension(dimension),centroid_tree(dimension){srand(time(0));}

    void initialize(vector<double*>& points) {
        centroids.clear();
        for (int i = 0; i < k; ++i) {
            centroids.push_back(points[rand() % points.size()]);
            centroid_tree.insert(centroids.back());
            cout<<centroids.back()[0]<<" - "<<centroids.back()[1]<<endl;
        }
    }

    double* closest_centroid(double*point){
        return centroid_tree.nearest(point);
    }

      void update_centroids(map<double*, vector<double*>>& clusters) {
        vector<double*>new_centroids;
        for (int i;i<k;i++){centroids.pop_back();}
        for (auto& cluster : clusters) {
            double* centroid = cluster.first;
            vector<double*>& points = cluster.second;
            if (points.empty()) continue;
            double* new_centroid = new double[dimension]();
            for (auto& point : points) {
                for (int d = 0; d < dimension; ++d) {
                    new_centroid[d] += point[d];
                }
            }
            for (int d = 0; d < dimension; ++d) {
                new_centroid[d] /= points.size();
            }
            centroids.push_back(new_centroid);
        }
        centroid_tree = KDtree(dimension);//reiniciamos kdtree
        for (auto& centroid : centroids) {
            centroid_tree.insert(centroid);
        }
    }



    void fit(vector<double*>& points){
        initialize(points);
        bool converged=false;
        map<double*, int> point_to_cluster;
        while(!converged){
            map<double*, vector<double*>> clusters;
            for(auto& point: points){
                double* centroid=closest_centroid(point);
                clusters[centroid].push_back(point);
                point_to_cluster[point] = centroid_to_cluster_id[centroid];
            }
            vector<double*> old_centroids=centroids;
            update_centroids(clusters);
            converged=true;
            for(int i=0;i<k;i++){
                if(centroid_tree.euclidian_dist(old_centroids[i],centroids[i])>0.0001){
                    converged=false;
                    break;
                }
            }
        }
        save_clusters_to_csv(point_to_cluster);
    }

    void save_clusters_to_csv(map<double*, int>& point_to_cluster) {
        ofstream file("clusters.csv");
        file << "x,y,cluster_id" << endl;
        for (const auto& entry : point_to_cluster) {
            double* point = entry.first;
            int cluster_id = entry.second;
            file << point[0] << "," << point[1] << "," << cluster_id << endl;
        }
        file.close();
    }

    void print_centroids(){
        centroid_tree.print_node(centroid_tree.root,0);
    }
};


int main() {
    int dimension = 2; 
    string filename = "data2k.csv"; 
    vector<double*> points = load_points_from_csv(filename, dimension);
    int k = 150; // Num de centroides 
    Kmeans kmeans(k, dimension);
    kmeans.fit(points);
    kmeans.print_centroids();
    

    for (auto& point : points) {
        delete[] point;
    }
    return 0;
}