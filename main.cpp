#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>
#include <algorithm>

std::vector<std::string> species = {"hg38", "panTro4", "gorGor3", "ponAbe2", "nomLeu3", "rheMac3", "macFas5", "papAnu2", "chlSab2", "calJac3", "saiBol1", "otoGar3", "tupChi1", "speTri2", "jacJac1", "micOch1", "criGri1", "mesAur1", "mm10", "rn6", "hetGla2", "cavPor3", "chiLan1", "octDeg1", "oryCun2", "ochPri3", "susScr3", "vicPac2", "camFer1", "turTru2", "orcOrc1", "panHod1", "bosTau8", "oviAri3", "capHir1", "equCab2", "cerSim1", "felCat8", "canFam3", "musFur1", "ailMel1", "odoRosDiv1", "lepWed1", "pteAle1", "pteVam1", "eptFus1", "myoDav1", "myoLuc2", "eriEur2", "sorAra2", "conCri1", "loxAfr3", "eleEdw1", "triMan1", "chrAsi1", "echTel2", "oryAfe1", "dasNov3", "monDom5", "sarHar1", "macEug2", "ornAna1", "colLiv1", "falChe1", "falPer1", "ficAlb2", "zonAlb1", "geoFor1", "taeGut2", "pseHum1", "melUnd1", "amaVit1", "araMac1", "anaPla1", "galGal4", "melGal1", "allMis1", "cheMyd1", "chrPic2", "pelSin1", "apaSpi1", "anoCar2", "xenTro7", "latCha1", "tetNig2", "fr3", "takFla1", "oreNil2", "neoBri1", "hapBur1", "mayZeb1", "punNye1", "oryLat2", "xipMac1", "gasAcu1", "gadMor1", "danRer10", "astMex1", "lepOcu1", "petMar2"};
std::string ref_species = "hg38";

class node {
private:
    node *parent;
    node *left_child;
    node *right_child;
    double branch_length;
    std::string name;

public:
    node(){
        this->branch_length = 0;
        this->name = "";
        parent = NULL;
        left_child = NULL;
        right_child = NULL;
    }

    void set_parent(node *parent) {
        this->parent = parent;
    }

    void set_left_child(node *left_child) {
        this->left_child = left_child;
    }

    void set_right_child(node *right_child) {
        this->right_child = right_child;
    }

    void set_name (std::string name) {
        this->name = name;
    }

    void set_branch_length(double bl) {
        this->branch_length = bl;
    }

    double get_branch_length() {
        return this->branch_length;
    }

    node* get_left_child() {
        return this->left_child;
    }

    node* get_right_child() {
        return this->right_child;
    }

    node* get_parent() {
        return this->parent;
    }

    std::string get_name() {
        return this->name;
    }
};

void _newick_parse(std::string& str, node* curr_node)
{
    if (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
    {
        while (str[0] != '(' && str[0] != ',' && str[0] != ')' && str[0] != ':')
        {
            curr_node->set_name(curr_node->get_name() += str[0]);
            str = str.substr(1); // TODO: inefficient!
        }
        if (find(species.begin(), species.end(), curr_node->get_name()) == species.end()) {
            curr_node->set_name("");
        }

        str = str.substr(1); // remove ':'

        std::string len = "";
        while (std::isdigit(str[0]) || str[0] == '.')
        {
            len += str[0];
            str = str.substr(1);
        }
        curr_node->set_branch_length(std::stod(len));
    }

    if (str[0] == '(')
    {
        node *left = new node;
        left->set_parent(curr_node);
        curr_node->set_left_child(left);
        //
        str = str.substr(1); // remove '(' // TODO: inefficient
        _newick_parse(str, left);

        node* right = new node;
        right->set_parent(curr_node);
        curr_node->set_right_child(right);
        _newick_parse(str, right);
    }

    if (str[0] == ',')
    {
        str = str.substr(1); // remove ','
        return; // go up one
    }

    if (str[0] == ')')
    {
        str = str.substr(1); // remove ')'

        if (str.size() == 0 || str[0] == ';')
            return;

        while (str[0] != ':') {
            str = str.substr(1);
        }
        str = str.substr(1); // remove ':'

        std::string len = "";
        while (std::isdigit(str[0]) || str[0] == '.')
        {
            len += str[0];
            str = str.substr(1);
        }
        curr_node->get_parent()->set_branch_length(std::stod(len));
    }
}

node* find_ref_species(std::string ref_species, node* cur_node) {
    if (cur_node->get_name() == ref_species) {
        return cur_node;
    }
    if (cur_node->get_left_child() != NULL) {
        node *left = find_ref_species(ref_species, cur_node->get_left_child());
        if (left != NULL) {
            return left;
        }
    }
    if (cur_node->get_right_child() != NULL) {
        node *right = find_ref_species(ref_species, cur_node->get_right_child());
        if (right != NULL) {
            return right;
        }
    }
    return NULL;
}

void change_topology(node *cur_node) {
    if (cur_node->get_parent() == NULL) {
        return;
    } else if (cur_node->get_parent()->get_left_child() != cur_node) {
        cur_node->get_parent()->set_right_child(cur_node->get_parent()->get_left_child());
        cur_node->get_parent()->set_left_child(cur_node);
    }
    change_topology(cur_node->get_parent());
}

void print_tree(std::ofstream &out, node* cur_node) {
    if (cur_node->get_left_child() != NULL) {
        out << "(";
        print_tree(out, cur_node->get_left_child());
    }
    if (cur_node->get_right_child() != NULL) {
        print_tree(out, cur_node->get_right_child());
        out << ")";
    }
    if (cur_node->get_parent() != NULL) {
        out << cur_node->get_name() << ":" << std::to_string(cur_node->get_branch_length());
        if (cur_node->get_parent()->get_left_child() == cur_node) {
            out << ",";
        }
    }
}

int main() {

    std::ifstream newick;
    newick.open("/Users/sukhwanpark/Downloads/100vertebrates.nh");
    std::string lines;
    std::string line;
    getline(newick, line);
    while (!line.empty()) {
        lines += line;
        while (lines.find(' ') != std::string::npos) {
            size_t pos = lines.find(' ');
            lines.erase(lines.begin() + pos);
        }
        getline(newick, line);
    }
    newick.close();
    node *first = new node;

    _newick_parse(lines, first);

    node *top = new node;
    node *top_right = new node;
    top->set_right_child(top_right);
    top_right->set_parent(top);
    bool top_is_root = false;
    if (lines != ";") {
        top->set_left_child(first);
        first->set_parent(top);
        std::string len;
        _newick_parse(lines, top_right);
        top_is_root = true;
    }

    // read is over, now change the tree and print
    node *ref_node;
    if (top_is_root == true) {
        ref_node = find_ref_species(ref_species, top);
    } else {
        ref_node = find_ref_species(ref_species, first);
    }

    change_topology(ref_node);

    //reroot the tree to ref_species
    /*node *new_top = new node;
    new_top->set_left_child(ref_node);
    new_top->set_right_child(ref_node->get_parent());
    double length = ref_node->get_branch_length();
    double leng_temp = ref_node->get_parent()->get_branch_length();
    ref_node->set_branch_length(length / 2);
    ref_node->get_parent()->set_branch_length(length / 2);
    node *parent = ref_node->get_parent();
    ref_node->set_parent(new_top);
    node *prev_node = new_top;
    while (parent->get_parent() != NULL || parent->get_parent()->get_parent() != NULL) {
        parent->set_left_child(parent->get_right_child());
        parent->set_right_child(parent->get_parent());
        length = parent->get_parent()->get_branch_length();
        parent->get_parent()->set_branch_length(leng_temp);
        leng_temp = length;
        parent = parent->get_parent();
        parent->get_left_child()->set_parent(prev_node);
        prev_node = parent->get_left_child();
    }
    parent->set_left_child(parent->get_right_child());
    */

    //now print the result
    std::ofstream newick_out;
    newick_out.open("/Users/sukhwanpark/Downloads/100vertebrates_test.nh");
    double sum = 0;
    if (top_is_root == true) {
        sum += top->get_left_child()->get_branch_length();
        sum += top->get_right_child()->get_branch_length();
        sum /= 2;
        top->get_left_child()->set_branch_length(sum);
        top->get_right_child()->set_branch_length(sum);
        print_tree(newick_out, top);
    } else {
        sum += first->get_left_child()->get_branch_length();
        sum += first->get_right_child()->get_branch_length();
        sum /= 2;
        first->get_left_child()->set_branch_length(sum);
        first->get_right_child()->set_branch_length(sum);
        print_tree(newick_out, first);
    }
    newick_out << ";";
    newick_out.close();

    return 0;
}
