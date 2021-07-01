#ifndef ATTRIBUTEDLOCALMOTIFCLUSTER_GRAPH_H
#define ATTRIBUTEDLOCALMOTIFCLUSTER_GRAPH_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <queue>
#include <fstream>
#include <set>
#include <math.h>

static int NibbleCompare(std::pair<int, double> t1, std::pair<int, double> t2){
    return t1.second > t2.second;
}

class Graph{
public:
    std::vector<int> qNode;
    std::vector<int> QueryAttributes; // note that the size of QueryAttributes cannot be too large

    //the edge list has to be preprocessed so that the NodeId are between 0 and (VNum-1)
    std::vector<std::vector<int>* > EdgeList;
    std::vector<std::vector<int>* > Attributes;

    std::vector<std::vector<double>* > OverlayAttrWeight;
    std::vector<double> OverlayAttrDegree;
    std::vector<std::pair<int, double>> OverlayNormalizedAttrPPR;
    unsigned int OverlayAttrNibbleIdx;
    
    std::vector<std::vector<int>* > label_prop;

    void Reset(){
        for(auto oaw : OverlayAttrWeight) delete(oaw);
        OverlayAttrWeight.clear();
        OverlayAttrDegree.clear();
        OverlayNormalizedAttrPPR.clear();
    }

    void QReset(){
        QueryAttributes.clear();
        qNode.clear();
    }

    // this->EdgeList and this->Attributes is initialized
    Graph(int VNum, std::string edgelist_aligned, std::string attributes_aligned){
        bool EdgeRepeat = false;
        this->EdgeList.reserve(VNum);
        this->Attributes.reserve(VNum);
        for(int i=0;i<VNum;i++){
            this->EdgeList.push_back(new std::vector<int>());
            this->Attributes.push_back(new std::vector<int>());
        }

        FILE * edgelist_in;
        std::ifstream attributes_in;

        edgelist_in = fopen(edgelist_aligned.c_str(), "r");
        attributes_in.open(attributes_aligned);

        int t1, t2;
        while (fscanf(edgelist_in, "%d%d", &t1, &t2) != EOF){
            if(t1 == t2) continue;

            this->EdgeList[t1]->push_back(t2);
            if (!EdgeRepeat)
                this->EdgeList[t2]->push_back(t1);
        }
        fclose(edgelist_in);

        int NodeId=0;
        std::string attributes_line;
        while (getline(attributes_in, attributes_line)){
            std::string::size_type temp_pos;
            int size = attributes_line.size();
            for(unsigned int i = 0; i<size;i++){
                temp_pos = attributes_line.find(' ',i);
                if(temp_pos<size){
                    int sub = stoi(attributes_line.substr(i,temp_pos-i));
                    if(Attributes[NodeId]->empty() || Attributes[NodeId]->back() < sub)
                        this->Attributes[NodeId]->push_back(sub);
                    i=temp_pos;
                }
            }
            NodeId++;
        }
        attributes_in.close();

        for(auto list : EdgeList)
            std::sort(list->begin(), list->end());
    }

    // this->QueryAttributes and this->qNode is initialized
    void SetQueryParameter(std::vector<int>* qAttribtues, std::vector<int>* QNode){
        this->QueryAttributes = std::vector<int>();
        int qAttributesNum = qAttribtues->size();
        for(unsigned int i=0;i<qAttributesNum;i++)
            QueryAttributes.push_back(qAttribtues->at(i));

        int qNodeNum = QNode->size();
        for(unsigned int i=0;i<qNodeNum;i++)
            qNode.push_back(QNode->at(i));
    }
    
    void LabelPropagationReverse(){ //newly added
        for(int i=0;i<QueryAttributes.size();i++){
            for(unsigned int j=0;j<label_prop[i]->size();j++){
                int prop_node = label_prop[i]->at(j);
                int prop_label = QueryAttributes[i];
                Attributes[prop_node]->erase(
                        std::remove(Attributes[prop_node]->begin(),Attributes[prop_node]->end(),prop_label),
                        Attributes[prop_node]->end());
            }
        }

        for(auto lpp: label_prop) delete(lpp);
        label_prop.clear();
    }

    void LabelPropagation(double over_density_thres, int query_node, int VNum, double epsilon){
    	double teleportation = 0.1;
        std::vector<double>  propfreqency;
        std::vector<std::vector<bool>* > query_labeled;

        for(int i=0;i<QueryAttributes.size();i++){
            label_prop.push_back(new std::vector<int>());
            propfreqency.push_back(0.0);
            query_labeled.push_back(new std::vector<bool>());
            for(int j=0;j<EdgeList.size();j++)
                query_labeled[query_labeled.size()-1]->push_back(false);
        }

        std::vector<double> ppr;
        ppr.reserve(VNum);

        //computing the ppr distribution
        std::vector<bool> in_push_nodes;
        in_push_nodes.reserve(VNum);

        std::vector<double> residual;
        residual.reserve(VNum);

        for(int i=0;i<VNum;i++){
            if(i==query_node){
                residual.push_back(1.0);
                in_push_nodes.push_back(true);
            }
            else {
                residual.push_back(0);
                in_push_nodes.push_back(false);
            }
        }

        std::vector<int> push_nodes; // nodes that is going to be local-pushed
        push_nodes.reserve(VNum);
        push_nodes.push_back(query_node);

        for(int i=0;i<VNum;i++)
            ppr.emplace_back(0.0);

        int size=0;
        while((size = push_nodes.size())>0){
            int pn = push_nodes[size-1];
            push_nodes.pop_back();
            in_push_nodes[pn] = false;

            double residual_pushed = 0.0;
            while(( residual[pn]/EdgeList[pn]->size()) > epsilon){
                ppr[pn] += teleportation * residual[pn];
                residual[pn] = (1-teleportation)*residual[pn] / 2.0;
                residual_pushed += residual[pn];
            }
            if(residual_pushed > 0){
                for(unsigned int i=0;i<EdgeList[pn]->size();i++) {
                    int neighbor = EdgeList[pn]->at(i);
                    residual[neighbor] += (residual_pushed / EdgeList[pn]->size());
                    if (!in_push_nodes[neighbor] && residual[neighbor] / EdgeList[neighbor]->size() > epsilon) {
                        push_nodes.push_back(neighbor);
                        in_push_nodes[neighbor] = true;
                    }
                }
            }
        }
        //computing the ppr distribution END

        for(unsigned int i=0;i<Attributes.size();i++)
            for(unsigned int j=0;j<Attributes[i]->size();j++){
                int label = Attributes[i]->at(j);
                for(int m=0;m<QueryAttributes.size();m++)
                    if(label == QueryAttributes[m]){
                        propfreqency[m] += ppr[i];
                        query_labeled[m]->at(i) = true;
                    }
            }

        for(int i=0;i<QueryAttributes.size();i++)
            propfreqency[i] *= over_density_thres;

        for(int qa=0;qa<QueryAttributes.size();qa++)
            for(unsigned int n=0;n<EdgeList.size();n++)
                if(!query_labeled[qa]->at(n)){
                    double labeled_nbr = 0.0;
                    for(unsigned int e=0;e<EdgeList[n]->size();e++){
                        int nbr = EdgeList[n]->at(e);
                        if(query_labeled[qa]->at(nbr))
                            labeled_nbr+= 1.0;
                    }
                    if(labeled_nbr/EdgeList[n]->size() > propfreqency[qa]){
                        label_prop[qa]->push_back(n);
                        Attributes[n]->push_back(QueryAttributes[qa]);
                        std::sort(Attributes[n]->begin(),Attributes[n]->end());
                    }
                }

        for(auto x: query_labeled) delete(x);  // prop version 1
    }

    void OverlayAttrSweep(int VNum){ // used in OverlayAttrNibble()
        for(int i=0;i<VNum;i++)
            if(EdgeList[i]->size() > 0) OverlayNormalizedAttrPPR[i].second /= EdgeList[i]->size();

        double vol=0 , total_vol=0, cut=0;
        double conductance=1.0;
        this->OverlayAttrNibbleIdx = 0;
        
        std::sort(OverlayNormalizedAttrPPR.begin(), OverlayNormalizedAttrPPR.end(), NibbleCompare);

        for(int i=0;i<VNum;i++){
            if(OverlayAttrDegree[OverlayNormalizedAttrPPR[i].first] > 0)
                total_vol += OverlayAttrDegree[OverlayNormalizedAttrPPR[i].first];
            else break;
        }

        std::vector<bool> InCluster;
        InCluster.reserve(VNum);
        for(int i=0;i<VNum;i++) InCluster.push_back(false);

        int SizeThreshold = (VNum > 1000)? 1000 : VNum;
        for(unsigned int i=0;i<SizeThreshold;i++){
            int node = OverlayNormalizedAttrPPR[i].first;
            if(OverlayAttrDegree[node] < 0) break;//////////////////######
            vol += OverlayAttrDegree[node];
            for(unsigned int j=0;j<EdgeList[node]->size();j++){
                if(InCluster[EdgeList[node]->at(j)])
                    cut -= OverlayAttrWeight[node]->at(j);
                else cut += OverlayAttrWeight[node]->at(j);
            }

            double temp_conductance;
            if(total_vol <= vol) temp_conductance = 1;
            else temp_conductance = (vol < (total_vol-vol))? cut/vol : cut/(total_vol-vol);

            if(temp_conductance < conductance){
                conductance = temp_conductance;
                this->OverlayAttrNibbleIdx = i;
            }
            InCluster[node] = true;
        }
    }

    void OverlayAttrNibble(int query_node, int VNum, double lambda, double epsilon){
        double teleportation = 0.1;
        this->OverlayAttrWeight.reserve(VNum);
        this->OverlayAttrDegree.reserve(VNum);
        for (unsigned int i = 0; i < VNum; i++){
            OverlayAttrWeight.push_back(new std::vector<double>());
            OverlayAttrWeight[i]->reserve(EdgeList[i]->size());
            for (unsigned int j = 0; j < EdgeList[i]->size(); j++)
                OverlayAttrWeight[i]->push_back(1-lambda);
            OverlayAttrDegree.push_back(-1);
        }
        std::vector<unsigned int> QueryMatchedAttributesNum;
        QueryMatchedAttributesNum.reserve(VNum);
        for(unsigned int i=0;i<VNum;i++){
            QueryMatchedAttributesNum.push_back(0);
            unsigned int count_q = 0;
            unsigned int count_attributes = 0;
            while(count_q < QueryAttributes.size() && count_attributes < Attributes[i]->size()){
                if(Attributes[i]->at(count_attributes) == QueryAttributes[count_q]){
                    QueryMatchedAttributesNum[i]++;
                    count_attributes++;
                    count_q++;
                }
                else if(Attributes[i]->at(count_attributes) < QueryAttributes[count_q])
                    count_attributes++;
                else count_q++;
            }
        }

        std::vector<bool> in_push_nodes;
        in_push_nodes.reserve(VNum);

        std::vector<double> residual;
        residual.reserve(VNum);

        for(int i=0;i<VNum;i++){
            if(i == query_node){
                residual.push_back(1.0);
                in_push_nodes.push_back(true);
            }
            else {
                residual.push_back(0);
                in_push_nodes.push_back(false);
            }
        }

        std::vector<int> push_nodes; // nodes that is going to be local-pushed
        push_nodes.reserve(VNum);
        for(int qn : qNode) push_nodes.push_back(qn);

        OverlayNormalizedAttrPPR.reserve(VNum);
        for(int i=0;i<VNum;i++)
            OverlayNormalizedAttrPPR.emplace_back(std::pair<int, double>(i, 0.0));

        int size=0;
        while((size = push_nodes.size())>0){
            int pn = push_nodes[size-1];
            push_nodes.pop_back();
            in_push_nodes[pn] = false;

            if(OverlayAttrDegree[pn] < 0){
                OverlayAttrDegree[pn] = (1-lambda)*EdgeList[pn]->size();
                if(QueryMatchedAttributesNum[pn] > 0){
                    for(unsigned int i=0;i<EdgeList[pn]->size();i++){
                        int n1 = EdgeList[pn]->at(i);
                        if(QueryMatchedAttributesNum[n1] > 0){
                            unsigned int pn_nbr_count = 0;
                            unsigned int n1_nbr_count = 0;
                            unsigned int TNum = 0;
                            while(pn_nbr_count < EdgeList[pn]->size() && n1_nbr_count < EdgeList[n1]->size()){
                                if(EdgeList[pn]->at(pn_nbr_count) == EdgeList[n1]->at(n1_nbr_count)){
                                    TNum += QueryMatchedAttributesNum[pn]*QueryMatchedAttributesNum[n1]
                                            *QueryMatchedAttributesNum[EdgeList[pn]->at(pn_nbr_count)];
                                    pn_nbr_count++;
                                    n1_nbr_count++;
                                }else if(EdgeList[pn]->at(pn_nbr_count) < EdgeList[n1]->at(n1_nbr_count))
                                    pn_nbr_count++;
                                else n1_nbr_count++;
                            }
                            OverlayAttrDegree[pn] += lambda*TNum;
                            OverlayAttrWeight[pn]->at(i) += lambda*TNum;
                        }
                    }
                }
            }

            double residual_pushed = 0.0;
            while(( residual[pn]/EdgeList[pn]->size()) > epsilon){
                OverlayNormalizedAttrPPR[pn].second += teleportation * residual[pn];
                residual[pn] = (1-teleportation)*residual[pn] / 2.0;
                residual_pushed += residual[pn];
            }
            if(residual_pushed > 0){
                for(unsigned int i=0;i<EdgeList[pn]->size();i++) {
                    int neighbor = EdgeList[pn]->at(i);
                    residual[neighbor] += (residual_pushed * OverlayAttrWeight[pn]->at(i)) / OverlayAttrDegree[pn];
                    if (!in_push_nodes[neighbor] && residual[neighbor] / EdgeList[neighbor]->size() > epsilon) {
                        push_nodes.push_back(neighbor);
                        in_push_nodes[neighbor] = true;
                    }
                }
            }
        }
        OverlayAttrSweep(VNum);
    }

    /*Validation functions*/
    void ValidOverlayAttrNibble_ATCLScore(bool if_prop,int VNum, double lambda, double epsilon){
        std::vector<int>* multinode_founded_cmty=nullptr;
        for(int querynode_count = 0;querynode_count < qNode.size();querynode_count++) {
            if(if_prop) LabelPropagation(1.0, qNode[querynode_count], VNum, epsilon);

            this->OverlayAttrNibble(qNode[querynode_count], VNum, lambda, epsilon);
            unsigned int NibbleIdx;
            NibbleIdx = this->OverlayAttrNibbleIdx;

            double attr_density = 0.0;
            std::vector<std::pair<int, double>> found_cmty;
            std::vector<int> node_label;
            //node_label[i] > 0 means node i connect to only nodes InCluster ;
            //node_label[i] < 0 means node i connect to nodes outside InCluster;
            //node_label[i] = 0 means node i has been deleted
            //(abs(node_label[i])-1) means the number of attributes i match
            for (unsigned int i = 0; i <= NibbleIdx; i++) {
                node_label.emplace_back(0);
                found_cmty.emplace_back(std::pair<int, double>(this->OverlayNormalizedAttrPPR[i].first, 0.0));
            }

            std::vector<bool> InCluster;
            InCluster.reserve(VNum);
            for (unsigned int i = 0; i < VNum; i++) InCluster.push_back(false);
            for (unsigned int i = 0; i <= NibbleIdx; i++)
                InCluster[this->OverlayNormalizedAttrPPR[i].first] = true;

            std::vector<int> QACount;
            for (unsigned int i=0;i<QueryAttributes.size();i++)
                QACount.emplace_back(0);

            for (unsigned int i = 0; i <= NibbleIdx; i++) {
                double out_R = 0;
                int node = this->OverlayNormalizedAttrPPR[i].first;
                for (unsigned int j = 0; j < EdgeList[node]->size(); j++)
                    if (!InCluster[EdgeList[node]->at(j)]) out_R += OverlayAttrWeight[node]->at(j);
                if (out_R <= 0) node_label[i] = 1;
                else node_label[i] = -1;

                int attribute_num = 1;
                for (int j=0;j<QueryAttributes.size();j++) {
                    int qattr = QueryAttributes[j];
                    if (std::find(Attributes[node]->begin(), Attributes[node]->end(), qattr) != Attributes[node]->end()) {
                        attribute_num++;
                        QACount[j]++;
                    }
                }
                node_label[i] *= attribute_num;
                //attr_density += (attribute_num-1);
            }

            for(int i=0;i<QACount.size();i++)
                attr_density += QACount[i]*QACount[i];
            attr_density /= (NibbleIdx + 1);
            int remained_node = (NibbleIdx + 1);

            bool MinAttrBetaIsGroundTruth = false;
            while (!MinAttrBetaIsGroundTruth) {
                double vol = 0, cut = 0;
                double conductance;

                for (int i = 0; i <= NibbleIdx; i++) {
                    if (node_label[i] != 0) {
                        int node = this->OverlayNormalizedAttrPPR[i].first;
                        vol += OverlayAttrDegree[node];
                        for (unsigned int j = 0; j < EdgeList[node]->size(); j++)
                            if (!InCluster[EdgeList[node]->at(j)]) cut += OverlayAttrWeight[node]->at(j);
                    }
                }
                conductance = cut / vol;

                for (int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] < 0) {
                        double closeness_R;
                        double in_R = 0, out_R = 0;
                        int node = this->OverlayNormalizedAttrPPR[i].first;
                        for (unsigned int j = 0; j < EdgeList[node]->size(); j++) {
                            if (!InCluster[EdgeList[node]->at(j)]) out_R += OverlayAttrWeight[node]->at(j);
                            else in_R += OverlayAttrWeight[node]->at(j);
                        }

                        closeness_R = in_R / out_R;

                        double beta_R = closeness_R * (conductance + 1) + conductance - 1;
                        beta_R /= ((cut / out_R) - conductance * (closeness_R + 1));
                        beta_R += 1.0;

                        found_cmty[i].second = beta_R;
                    }

                for (int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] > 0) {
                        double temp_cut = cut, temp_vol = vol;
                        int node = this->OverlayNormalizedAttrPPR[i].first;
                        temp_cut += OverlayAttrDegree[node];
                        temp_vol -= OverlayAttrDegree[node];
                        double temp_conductance = temp_cut / temp_vol;
                        found_cmty[i].second = (temp_conductance / conductance);
                    }

                for (int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] != 0) {
                        //found_cmty[i].second *= (abs(node_label[i]) - 1);
                        int atc_label_score = 0;
                        for(int j=0;j<QueryAttributes.size();j++){
                            int ql = QueryAttributes[j];
                            int node = OverlayNormalizedAttrPPR[i].first;
                            if(std::find(Attributes[node]->begin(), Attributes[node]->end(), ql) != Attributes[node]->end()){
                                atc_label_score += 2*QACount[j]-1;
                            }
                        }
                        found_cmty[i].second *= atc_label_score;
                    }

                unsigned int min_beta_pos = 0;
                for (unsigned int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] != 0) {
                        min_beta_pos = i;
                        break;
                    }

                for (unsigned int i = 0; i <= NibbleIdx; i++) {
                    if (node_label[i] != 0 && (found_cmty[i].second - found_cmty[min_beta_pos].second) <= 0)
                        min_beta_pos = i;
                }

                int min_beta_node = found_cmty[min_beta_pos].first;

                double pre_attr_density = attr_density;
                //attr_density *= remained_node;

                if(remained_node <= 1)
                //if (remained_node == 1 || (abs(node_label[min_beta_pos]) - 1) >= this->QueryAttributes.size())
                    MinAttrBetaIsGroundTruth = true;
                else {
                    attr_density = 0.0;
                    for(int j=0;j<QACount.size();j++){
                        int ql = QueryAttributes[j];
                        if(std::find(Attributes[min_beta_node]->begin(), Attributes[min_beta_node]->end(), ql) != Attributes[min_beta_node]->end()){
                            QACount[j]--;
                        }
                        attr_density += QACount[j]*QACount[j];
                    }

                    //attr_density -= (abs(node_label[min_beta_pos]) - 1);
                    attr_density /= (remained_node - 1);
                    MinAttrBetaIsGroundTruth = (attr_density <= pre_attr_density);
                }

                if (!MinAttrBetaIsGroundTruth) {
                    node_label[min_beta_pos] = 0;
                    InCluster[min_beta_node] = false;
                    remained_node--;

                    for (unsigned int i = 0; i <= NibbleIdx; i++) {
                        if (node_label[i] != 0) {
                            double out_R = 0;
                            int node = this->OverlayNormalizedAttrPPR[i].first;
                            for (unsigned int j = 0; j < EdgeList[node]->size(); j++)
                                if (!InCluster[EdgeList[node]->at(j)]) out_R += OverlayAttrWeight[node]->at(j);
                            if (out_R <= 0) node_label[i] = abs(node_label[i]);
                            else node_label[i] = -abs(node_label[i]);
                        }
                    }
                }
            }

            auto current_founded_cmty=new std::vector<int>();
            std::vector<int>* intersected_cmty=nullptr;
            for(int i=0;i<=NibbleIdx;i++)
                if(node_label[i] != 0) current_founded_cmty->push_back(found_cmty[i].first);
            std::sort(current_founded_cmty->begin(), current_founded_cmty->end());

            if(querynode_count == 0) intersected_cmty = current_founded_cmty;
            else{
                intersected_cmty = new std::vector<int>();
                for(unsigned int multi_count = 0, current_count = 0;
                    multi_count < multinode_founded_cmty->size() && current_count<current_founded_cmty->size();){
                    if(multinode_founded_cmty->at(multi_count) == current_founded_cmty->at(current_count))
                    {
                        intersected_cmty->push_back(multinode_founded_cmty->at(multi_count));
                        multi_count++;current_count++;
                    }
                    else if(multinode_founded_cmty->at(multi_count) < current_founded_cmty->at(current_count))
                        multi_count++;
                    else current_count++;
                }
            }

            if(querynode_count != 0) delete(current_founded_cmty);
            delete(multinode_founded_cmty);
            multinode_founded_cmty = intersected_cmty;
            Reset();
            if(if_prop) LabelPropagationReverse();
        }

        for (unsigned int i = 0; i < multinode_founded_cmty->size(); i++) {
            std::cout<<multinode_founded_cmty->at(i)<<" ";
        }
        std::cout<<std::endl;
    }
    
    void ValidOverlayAttrNibble_SimpleLScore(bool if_prop,int VNum, double lambda, double epsilon){
        std::vector<int>* multinode_founded_cmty=nullptr;
        for(int querynode_count = 0;querynode_count < qNode.size();querynode_count++) {
            if(if_prop) LabelPropagation(1.0, qNode[querynode_count], VNum, epsilon);
            this->OverlayAttrNibble(qNode[querynode_count], VNum, lambda, epsilon);
            unsigned int NibbleIdx;
            NibbleIdx = this->OverlayAttrNibbleIdx;

            double attr_density = 0.0;
            std::vector<std::pair<int, double>> found_cmty;
            std::vector<int> node_label;
            //node_label[i] > 0 means node i connect to only nodes InCluster ;
            //node_label[i] < 0 means node i connect to nodes outside InCluster;
            //node_label[i] = 0 means node i has been deleted
            //(abs(node_label[i])-1) means the number of attributes i match
            for (unsigned int i = 0; i <= NibbleIdx; i++) {
                node_label.emplace_back(0);
                found_cmty.emplace_back(std::pair<int, double>(this->OverlayNormalizedAttrPPR[i].first, 0.0));
            }

            std::vector<bool> InCluster;
            InCluster.reserve(VNum);
            for (unsigned int i = 0; i < VNum; i++) InCluster.push_back(false);
            for (unsigned int i = 0; i <= NibbleIdx; i++)
                InCluster[this->OverlayNormalizedAttrPPR[i].first] = true;

            for (unsigned int i = 0; i <= NibbleIdx; i++) {
                double out_R = 0;
                int node = this->OverlayNormalizedAttrPPR[i].first;
                for (unsigned int j = 0; j < EdgeList[node]->size(); j++)
                    if (!InCluster[EdgeList[node]->at(j)]) out_R += OverlayAttrWeight[node]->at(j);
                if (out_R <= 0) node_label[i] = 1;
                else node_label[i] = -1;

                int attribute_num = 1;
                for (int qattr : QueryAttributes)
                    if (std::find(Attributes[node]->begin(), Attributes[node]->end(), qattr) != Attributes[node]->end())
                        attribute_num++;
                node_label[i] *= attribute_num;
                attr_density += (attribute_num - 1);
            }
            attr_density /= (NibbleIdx + 1);
            int remained_node = (NibbleIdx + 1);

            bool MinAttrBetaIsGroundTruth = false;
            while (!MinAttrBetaIsGroundTruth) {
                double vol = 0, cut = 0;
                double conductance;

                for (int i = 0; i <= NibbleIdx; i++) {
                    if (node_label[i] != 0) {
                        int node = this->OverlayNormalizedAttrPPR[i].first;
                        vol += OverlayAttrDegree[node];
                        for (unsigned int j = 0; j < EdgeList[node]->size(); j++)
                            if (!InCluster[EdgeList[node]->at(j)]) cut += OverlayAttrWeight[node]->at(j);
                    }
                }
                conductance = cut / vol;

                for (int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] < 0) {
                        double closeness_R;
                        double in_R = 0, out_R = 0;
                        int node = this->OverlayNormalizedAttrPPR[i].first;
                        for (unsigned int j = 0; j < EdgeList[node]->size(); j++) {
                            if (!InCluster[EdgeList[node]->at(j)]) out_R += OverlayAttrWeight[node]->at(j);
                            else in_R += OverlayAttrWeight[node]->at(j);
                        }

                        closeness_R = in_R / out_R;

                        double beta_R = closeness_R * (conductance + 1) + conductance - 1;
                        beta_R /= ((cut / out_R) - conductance * (closeness_R + 1));
                        beta_R += 1.0;

                        found_cmty[i].second = beta_R;
                    }

                for (int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] > 0) {
                        double temp_cut = cut, temp_vol = vol;
                        int node = this->OverlayNormalizedAttrPPR[i].first;
                        temp_cut += OverlayAttrDegree[node];
                        temp_vol -= OverlayAttrDegree[node];
                        double temp_conductance = temp_cut / temp_vol;
                        found_cmty[i].second = (temp_conductance / conductance);
                    }

                for (int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] != 0) {
                        found_cmty[i].second *= (abs(node_label[i]) - 1);
                    }

                unsigned int min_beta_pos = 0;
                for (unsigned int i = 0; i <= NibbleIdx; i++)
                    if (node_label[i] != 0) {
                        min_beta_pos = i;
                        break;
                    }

                for (unsigned int i = 0; i <= NibbleIdx; i++) {
                    if (node_label[i] != 0 && (found_cmty[i].second - found_cmty[min_beta_pos].second) <= 0)
                        min_beta_pos = i;
                }

                int min_beta_node = found_cmty[min_beta_pos].first;

                double pre_attr_density = attr_density;
                attr_density *= remained_node;

                if (remained_node == 1 || (abs(node_label[min_beta_pos]) - 1) >= this->QueryAttributes.size())
                    MinAttrBetaIsGroundTruth = true;
                else {
                    attr_density -= (abs(node_label[min_beta_pos]) - 1);
                    attr_density /= (remained_node - 1);
                    MinAttrBetaIsGroundTruth = (attr_density <= pre_attr_density);
                }

                if (!MinAttrBetaIsGroundTruth) {
                    node_label[min_beta_pos] = 0;
                    InCluster[min_beta_node] = false;
                    remained_node--;

                    for (unsigned int i = 0; i <= NibbleIdx; i++) {
                        if (node_label[i] != 0) {
                            double out_R = 0;
                            int node = this->OverlayNormalizedAttrPPR[i].first;
                            for (unsigned int j = 0; j < EdgeList[node]->size(); j++)
                                if (!InCluster[EdgeList[node]->at(j)]) out_R += OverlayAttrWeight[node]->at(j);
                            if (out_R <= 0) node_label[i] = abs(node_label[i]);
                            else node_label[i] = -abs(node_label[i]);
                        }
                    }
                }
            }

            auto current_founded_cmty=new std::vector<int>();
            std::vector<int>* intersected_cmty=nullptr;
            for(int i=0;i<=NibbleIdx;i++)
                if(node_label[i] != 0) current_founded_cmty->push_back(found_cmty[i].first);
            std::sort(current_founded_cmty->begin(), current_founded_cmty->end());

            if(querynode_count == 0) intersected_cmty = current_founded_cmty;
            else{
                intersected_cmty = new std::vector<int>();
                for(unsigned int multi_count = 0, current_count = 0;
                    multi_count < multinode_founded_cmty->size() && current_count<current_founded_cmty->size();){
                    if(multinode_founded_cmty->at(multi_count) == current_founded_cmty->at(current_count))
                    {
                        intersected_cmty->push_back(multinode_founded_cmty->at(multi_count));
                        multi_count++;current_count++;
                    }
                    else if(multinode_founded_cmty->at(multi_count) < current_founded_cmty->at(current_count))
                        multi_count++;
                    else current_count++;
                }
            }

            if(querynode_count != 0) delete(current_founded_cmty);
            delete(multinode_founded_cmty);
            multinode_founded_cmty = intersected_cmty;
            Reset();
            if(if_prop) LabelPropagationReverse();
        }

        for (unsigned int i = 0; i < multinode_founded_cmty->size(); i++) {
            std::cout<<multinode_founded_cmty->at(i)<<" ";
        }
        std::cout<<std::endl;
    }
    /*Validation functions ends*/
};

#endif