#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <set>
#include <climits>
#include <queue>
#include <stack>
using namespace std;

#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <fstream>
#include <memory>
#include <sstream>
#include <set>
#include <iostream>

unordered_map<int,set<int>> Components;

class CDR3Entry {
private:
    int id;                 
    std::string sequence;   
    std::string Vgene;      
    std::string Jgene; 
    int count;

public:
    CDR3Entry(int _id, const std::string& _sequence, 
              const std::string& _Vgene, const std::string& _Jgene, 
              int _count)
        : id(_id), sequence(_sequence), Vgene(_Vgene), Jgene(_Jgene), count(_count) {}
    
    int getId() const { return id; }
    const std::string& getSequence() const { return sequence; }
    const std::string& getVgene() const { return Vgene; }
    const std::string& getJgene() const { return Jgene; }
    int getCount() const { return count; }
    
    void Print_Entry() const {
        cout << "CDR3:" << id << " | V:" << Vgene << " | J:" << Jgene 
           << " | COUNT:" << count << " | SEQ:" << sequence << endl;
    }
};


struct VJPair {
    std::string Vgene;
    std::string Jgene;
    
    bool operator<(const VJPair& other) const {
        if (Vgene != other.Vgene)
            return Vgene < other.Vgene;
        return Jgene < other.Jgene;
    }
    
    VJPair(const std::string& v, const std::string& j) {
        Vgene = v;
        Jgene = j;
    }
    
    void print_VJ() const {
        cout << Vgene << ":" << Jgene << endl;
    }
};

class CDR3Dataset {
private:
    std::vector<std::shared_ptr<CDR3Entry>> entries;
    
    std::unordered_map<int, std::shared_ptr<CDR3Entry>> idMap;
    std::map<VJPair, std::vector<int>> VJMap;
    bool parseHeader(const std::string& header, int& id, std::string& Vgene, 
                    std::string& Jgene, int& count) {
        //>CDR3:1|V_hit:IGHV3-21*01|J_hit:IGHJ3*01|COUNT:69
        size_t idPos = header.find("CDR3:");
        size_t vPos = header.find("V_hit:");
        size_t jPos = header.find("J_hit:");
        size_t countPos = header.find("COUNT:");
        
        
        // ID
        size_t idEnd = header.find("|", idPos);
        if (idEnd == std::string::npos) idEnd = header.length();
        std::string idStr = header.substr(idPos + 5, idEnd - (idPos + 5));
        id = std::stoi(idStr);
        
        // V gene 
        size_t vEnd = header.find("*", vPos);
        if (vEnd == std::string::npos) vEnd = header.length();
        Vgene = header.substr(vPos + 6, vEnd - (vPos + 6));
        
        // J gene
        size_t jEnd = header.find("*", jPos);
        if (jEnd == std::string::npos) jEnd = header.length();
        Jgene = header.substr(jPos + 6, jEnd - (jPos + 6));
        
        // Count
        std::string countStr = header.substr(countPos + 6);
        count = std::stoi(countStr);
        
        return true;
    }
    
public:

    bool Read_Fasta(const std::string& filename) {
        std::ifstream file(filename);
        
        std::string line, header, sequence;
        int id = 0;
        std::string Vgene, Jgene;
        int count = 0;
        
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            
            if (line[0] == '>') {
                // add new entry
                if (!header.empty() && !sequence.empty()) {
                    if (parseHeader(header, id, Vgene, Jgene, count)) {
                        addEntry(id, sequence, Vgene, Jgene, count);
                    }
                }
                
                // collecting
                header = line;
                sequence.clear();
            } else {
                sequence += line;
            }
        }

        if (!header.empty() && !sequence.empty()) {
            if (parseHeader(header, id, Vgene, Jgene, count)) {
                addEntry(id, sequence, Vgene, Jgene, count);
            }
        }
        
        file.close();
        return true;
    }
    // function to add entry
    void addEntry(int id, const std::string& sequence, const std::string& Vgene, 
                  const std::string& Jgene, int count) {
        auto entry = std::make_shared<CDR3Entry>(id, sequence, Vgene, Jgene, count);
        
        entries.push_back(entry);
        idMap[id] = entry;
        
        VJPair vjPair(Vgene, Jgene);
        VJMap[vjPair].push_back(id);
    }
    
    std::shared_ptr<CDR3Entry> findById(int id) const {
        auto it = idMap.find(id);
        if (it != idMap.end()) {
            return it->second;
        }
        return nullptr;
    }
    
    std::vector<int> getCDR3IdsByVJPair(const std::string& Vgene, const std::string& Jgene) const {
        VJPair vjPair(Vgene, Jgene);
        auto it = VJMap.find(vjPair);
        if (it != VJMap.end()) {
            return it->second;
        }
        return {};
    }

    std::vector<VJPair> getAllVJPairs() const {
        std::vector<VJPair> result;
        for (const auto& pair : VJMap) {
            result.push_back(pair.first);
        }
        return result;
    }
    
};

int hammingDist(string str1, string str2) 
{ 
    int i = 0, count = 0; 
    while (str1[i] != '\0') { 
        if (str1[i] != str2[i]) 
            count++; 
        i++; 
    } 
    return count; 
} 

unordered_map<int,set<int>> Hamming_graph(const CDR3Dataset& data){
    unordered_map<int,set<int>> Vertex_map;
    std::vector<VJPair> ALLVJPairs = data.getAllVJPairs();
    for (const auto& vjpair : ALLVJPairs){
        vector<int> same_VJ_vec = data.getCDR3IdsByVJPair(vjpair.Vgene,vjpair.Jgene);
        for (int i=0; i < same_VJ_vec.size();i++){
            for (int j=i+1; j < same_VJ_vec.size();j++){
                // edge ruls
                int v = same_VJ_vec[i];
                int u = same_VJ_vec[j];
                string seq_v = data.findById(v)->getSequence();
                string seq_u = data.findById(u)->getSequence();
                if (seq_u.size() == seq_v.size()){
                    if ((double)hammingDist(seq_v,seq_u)/(double)seq_v.size() <= 0.1){
                        Vertex_map[v].insert(u);
                        Vertex_map[u].insert(v);
                    }
                }
            }
        }
    }
    return Vertex_map;
}

map<int,int> clonal_lineage_search(const unordered_map<int,set<int>>& Vertex_map){
    map<int,int> clonal_lineage;
    int component = 1;
    for (const auto& v: Vertex_map){
        if (clonal_lineage[v.first] == 0){
            queue<int> q;
            q.push(v.first);
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                Components[component].insert(u);
                
                for (int i : Vertex_map.at(u)) {
                    if (clonal_lineage[i] == 0) {
                        clonal_lineage[i] = component;
                        q.push(i);
                        //cout << i << " : " << component << endl;
                    }
                }
            }
            component++;
        }
    }
    return clonal_lineage;
}


int main() {
    CDR3Dataset data;
    data.Read_Fasta("./Alignment_results/compressed_cdr3s.fasta");
    unordered_map<int,set<int>> Vertex_map = Hamming_graph(data);
    // int num_edge = 0;
    // int num_vertex = 0;
    // for (const auto& v : Vertex_map){
    //     num_vertex++;
    //     for (auto u : v.second){
    //         num_edge++;
    //     }
    // }
    // cout << num_edge << " ; " << num_vertex << endl;
    map<int,int> clonal_lineage = clonal_lineage_search(Vertex_map);
    ofstream outFile("ground_truth.csv");
    cout << "Sequence_ID, Cluster" << endl;
    for (auto pair: clonal_lineage){
        outFile << pair.first << ", " << pair.second << endl;
    }
    outFile.close();
    return 0;
}