#include <time.h>
#include <unistd.h>
#include "Graph.h"


void SimpleLScore(std::string filename ,int VNum, double epsilon, double lambda, bool if_prop){
    std::string edgelist_aligned = "./datasets/";
    edgelist_aligned += filename;
    edgelist_aligned += ".edges";

    std::string attributes_aligned = "./datasets/";
    attributes_aligned += filename;
    attributes_aligned += ".feat";

    std::string query_aligned = "./datasets/";
    query_aligned += filename;
    query_aligned += ".query";

    auto g = new Graph(VNum, edgelist_aligned, attributes_aligned);

    int g_qnum = 0;
    std::ifstream query_in;
    query_in.open(query_aligned);
    std::string query_line;

    while(getline(query_in,query_line)){
        auto QNode = new std::vector<int>();

        int QNodeCount=0;
        std::string::size_type temp_pos;
        for(unsigned int i=0;i<query_line.size();i++){
            temp_pos = query_line.find(' ',i);
            if(temp_pos < query_line.size()){
                int sub = stoi(query_line.substr(i, temp_pos-i));
                QNode->push_back(sub);
                i = temp_pos;
            }
            QNodeCount++;
        }

        getline(query_in,query_line);

        auto qAttributes = new std::vector<int>();
        for(unsigned int i=0;i<query_line.size();i++){
            temp_pos = query_line.find(' ',i);
            if(temp_pos < query_line.size()){
                int sub = stoi(query_line.substr(i, temp_pos-i));
                qAttributes->push_back(sub);
                i = temp_pos;
            }
            QNodeCount++;
        }
        std::sort(qAttributes->begin(), qAttributes->end());

        std::cout<<"#"<<g_qnum<<" ";

        g->SetQueryParameter(qAttributes, QNode);
        g->ValidOverlayAttrNibble_SimpleLScore(false,VNum, lambda, epsilon);
        g->QReset();
        g_qnum++;
        delete(qAttributes);
        delete(QNode);
    }
    query_in.close();
}

void ATCLScore(std::string filename ,int VNum, double epsilon, double lambda, bool if_prop) {
    std::string edgelist_aligned = "./datasets/";
    edgelist_aligned += filename;
    edgelist_aligned += ".edges";

    std::string attributes_aligned = "./datasets/";
    attributes_aligned += filename;
    attributes_aligned += ".feat";

    std::string query_aligned = "./datasets/";
    query_aligned += filename;
    query_aligned += ".query";

    auto g = new Graph(VNum, edgelist_aligned, attributes_aligned);

    int g_qnum = 0;
    std::ifstream query_in;
    query_in.open(query_aligned);
    std::string query_line;

    while(getline(query_in,query_line)){
        auto QNode = new std::vector<int>();

        int QNodeCount=0;
        std::string::size_type temp_pos;
        for(unsigned int i=0;i<query_line.size();i++){
            temp_pos = query_line.find(' ',i);
            if(temp_pos < query_line.size()){
                int sub = stoi(query_line.substr(i, temp_pos-i));
                QNode->push_back(sub);
                i = temp_pos;
            }
            QNodeCount++;
        }

        getline(query_in,query_line);

        auto qAttributes = new std::vector<int>();
        for(unsigned int i=0;i<query_line.size();i++){
            temp_pos = query_line.find(' ',i);
            if(temp_pos < query_line.size()){
                int sub = stoi(query_line.substr(i, temp_pos-i));
                qAttributes->push_back(sub);
                i = temp_pos;
            }
            QNodeCount++;
        }
        std::sort(qAttributes->begin(), qAttributes->end());

        std::cout<<"#"<<g_qnum<<" ";

        g->SetQueryParameter(qAttributes, QNode);
        g->ValidOverlayAttrNibble_ATCLScore(false,VNum, lambda, epsilon);
        g->QReset();
        g_qnum++;
        delete(qAttributes);
        delete(QNode);
    }
    query_in.close();
}

int main(int argc, char **argv){
    int opt;
    const char *optstring = "f:n:l:e:p:s:";
    //-f: name of the dataset
    //-n: number of nodes in the network
    //-l: lambda
    //-e: error tolerance epsilon for PPR
    //-p: if label propagation preprocessing is applied
    //-s: whether the atc label score function or the simple label score function is used
    
    std::string filename;
    int VNum;
    double epsilon;
    double lambda;
    bool if_prop;
    int score_function_type;
    while ((opt=getopt(argc, argv, optstring)) != -1){
        switch(opt){
            case 'f':
                filename = optarg;
                break;
            case 'n':
                VNum = atoi(optarg);
                break;
            case 'l':
                lambda = atof(optarg);
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'p':
                if_prop = atoi(optarg);
                break;
            case 's':
                score_function_type = atoi(optarg);
                break;
        }
    }
    
    switch(score_function_type){
    	case 0:
    		SimpleLScore(filename ,VNum, epsilon, lambda, if_prop);break;
     	case 1:
    		ATCLScore(filename ,VNum, epsilon, lambda, if_prop);break;
    }
    return 0;
}