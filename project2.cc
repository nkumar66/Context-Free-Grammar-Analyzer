/*
 * Copyright (C) Mohsen Zohrevandi, 2017
 *               Rida Bazzi 2019
 * Do not share this file with anyone
 */

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "lexer.h"

using namespace std;

struct GrammarRule {
    string lhs;
    vector<string> rhs;
};

LexicalAnalyzer lexer;

vector<GrammarRule> grammar_rules;
vector<string> terminals;
vector<string> non_terminals;
vector<string> unconfirmed; //for if we don't know yet, to ensure that they are put in the correct order

unordered_set<string> seen_terminals;   //this needs to be an unordered set to ensure unique insertions
unordered_set<string> seen_non_terminals;   //this needs to be an unordered set to ensure unique insertions
unordered_set<string> nullable_set;
string current_lhs;
unordered_map<string, set<string>> first_sets;  //for task 3
unordered_map<string, set<string>> follow_sets; //for task 4



void Grammar();
void Rule_list();
void Rule();
void Right_hand_side();

void syntax_error() {
    cout << "SYNTAX ERROR !!!!!!!!!!!!!!\n";
    exit(1);
}

Token expect(TokenType expected_type) {
    Token t = lexer.GetToken();
    if (t.token_type != expected_type) syntax_error();
    return t;
}

//function to see if an element already exists in the vector
bool contains(const vector<string>& vec, const string& str) {
    return find(vec.begin(), vec.end(), str) != vec.end();

}

void add_non_terminal(const string& s) {
    //if the terminal is not found in the set already, then add it in order
    if (!contains(non_terminals, s)) {
        non_terminals.push_back(s);
        seen_non_terminals.insert(s);

        unconfirmed.push_back(s);
    }
}

void add_terminal(const string& s) {
    //if the terminal is not found in the set already, then add it in order
    if (!contains(terminals, s)) {
        terminals.push_back(s);
        seen_terminals.insert(s);

        unconfirmed.push_back(s);
    }
}



void add_unconfirmed(const string& name) {

    //if the terminal is not found in the set already, then add it in order
    if (!contains(unconfirmed, name))
        unconfirmed.push_back(name);
}

//because of when i have wrong terminals
void fix_terminals() {
    //go through all the symbols in terminals, and check if it already exists in seen non terminals, and if it doesn't exist there, add it to fixed_terminals since it'a  terminal
    vector<string> fixed;
    for (const auto& t : terminals) {
        if (seen_non_terminals.find(t) == seen_non_terminals.end())
            fixed.push_back(t);
    }
    terminals = fixed;

    //now we have to remove elements from unconfirmed that appear in terminals
    vector<string> cleaned;
    unordered_set<string> added;
    for (const auto& s : unconfirmed) {
        //only keep symbols that arne't in the terminals
        if (!contains(terminals, s) && added.find(s) == added.end()) {
            cleaned.push_back(s);
            added.insert(s);

        }
    }
    unconfirmed = cleaned;
    non_terminals = unconfirmed;

}





void ReadGrammar() {
    //call the Grammar function, to start the parsing processs
    Grammar();
}


void Grammar() {
    //Grammar → Rule-list HASH
    Token t = lexer.peek(1);
    if (t.token_type == ID) {
        Rule_list();
        expect(HASH);
    } else {
        syntax_error();
    }
}

void Rule_list() {
    //Rule-list → Rule Rule-list | Rule
    Rule();
    Token t = lexer.peek(1);
    if (t.token_type == ID) {
        Rule_list();
    }

}

//modified to add a parameter because it wasn't working earlier when it took no parameters, just so that I can add
//everything into the grammar rules correctly, and it doesn't skip over anything, make sit easier
void Id_list(vector<string>& rhs) {
    //Id-list → ID Id-list | e
    Token t = lexer.peek(1);
    if (t.token_type == ID) {
        
        t = expect(ID);
        rhs.push_back(t.lexeme);

        
        if (seen_non_terminals.find(t.lexeme) != seen_non_terminals.end()) {
            add_unconfirmed(t.lexeme);
        } else {
            add_terminal(t.lexeme);
            add_unconfirmed(t.lexeme);
        }
        Id_list(rhs);

    }
}

void Rule() {
    //Rule → ID ARROW Right-hand-side STAR
    Token lhs = expect(ID);
    add_non_terminal(lhs.lexeme);

    expect(ARROW);
    current_lhs = lhs.lexeme;   //to add the lhs to store it, we will need it to reference later

    Right_hand_side();
    expect(STAR);
}


void Right_hand_side() {
    //Right-hand-side → Id-list | Id-list OR Right-hand-side
    vector<string> rhs;
    Token t = lexer.peek(1);

    if (t.token_type == ID) {
        //if its not a non terminal, then it must be a terminal so lets check if it's in the non temrinals, and if not, then add it to the terminals
        //rhs.push_back(t.lexeme);
        if (seen_non_terminals.find(t.lexeme) == seen_non_terminals.end()) {
            add_terminal(t.lexeme);
            add_unconfirmed(t.lexeme);
        } else {
            add_unconfirmed(t.lexeme);
        }
        Id_list(rhs);
    }

    grammar_rules.push_back({current_lhs, rhs});

    t = lexer.peek(1);
    if (t.token_type == OR) {
        expect(OR);
        Right_hand_side();
    }
}



void compute_nullable() {

    //we need to go through the grammar rules and if the right hand side is empty, add the left hand side
    for (const auto& rule : grammar_rules) {
        if (rule.rhs.empty()) {
            nullable_set.insert(rule.lhs);
        }
    }

    bool changed;
    do {
        changed = false;
        //go through teh grammar rules and skip if it already exists
        for (const auto& rule : grammar_rules) {
            if (nullable_set.find(rule.lhs) != nullable_set.end()) continue;

            bool all_nullable = true;
            //now we have to check if all the rules in the right hand side are nullable or not
            for (const string& s : rule.rhs) {
                if (seen_non_terminals.find(s) == seen_non_terminals.end() || nullable_set.find(s) == nullable_set.end()) {
                    all_nullable = false;
                    break;
                }
            }

            if (all_nullable) {
                nullable_set.insert(rule.lhs);
                changed = true;
            }
        }
    } while (changed);
}

/* 
 * Task 1: 
 * Printing the terminals, then nonterminals of grammar in appearing order
 * output is one line, and all names are space delineated
*/
void Task1() {
    fix_terminals();
    for (const auto& t : terminals) cout << t << " ";
    for (const auto& nt : non_terminals) cout << nt << " ";
    cout << endl;
}


/*
 * Task 2:
 * Print out nullable set of the grammar in specified format.
*/
void Task2() {
    compute_nullable();
    cout << "Nullable = {";
    bool first = true;
    set<string> printed;

    for (const string& symbol : unconfirmed) {
        if (seen_non_terminals.count(symbol) && nullable_set.count(symbol) && !printed.count(symbol)) {
            if (!first) cout << ",";
            cout << symbol;
            printed.insert(symbol);
            first = false;
        }
    }

    cout << "}" << endl;
}

// Task 3: FIRST sets

void compute_FIRST() {
    // before anything, lets ensure that the classification of symbols is fixed
    fix_terminals();

    // initialize FIRST sets
    for (const string & symbol : terminals) {
        first_sets[symbol] = { symbol };
    }
    // for non-terminals, start with an empty FIRST set
    for (const string & nt : non_terminals) {
        first_sets[nt] = {};
    }

    bool changed;
    do {
        changed = false;
        for (const auto & rule : grammar_rules) {
            string A = rule.lhs;
            //process the RHS one symbol at a time from rule.rhs
            for (const string & X : rule.rhs) {
                //if X is a terminal, add it to FIRST(A) and stop for this production
                if (find(terminals.begin(), terminals.end(), X) != terminals.end()) {
                    if (first_sets[A].find(X) == first_sets[A].end()) {
                        first_sets[A].insert(X);

                        changed = true;
                    }
                    // terminals are never going to be n ullable, so we stop here
                    break;
                }
                // else X is a non-terminal.
                if (find(non_terminals.begin(), non_terminals.end(), X) != non_terminals.end()) {
                    // add all symbols from FIRST(X) to FIRST(A)
                    for (const string & t : first_sets[X]) {
                        if (first_sets[A].find(t) == first_sets[A].end()) {
                            first_sets[A].insert(t);
                            changed = true;
                        }
                    }
                    // If X is not nullable, then stop processing further symbols.
                    if (nullable_set.find(X) == nullable_set.end()) {
                        break;
                    }
                    // Otherwise, if X is nullable, continue to the next symbol.
                }
            }

        }
    } while (changed);
}

void Task3() {
    // ix terminals and non-terminals before computing FIRST sets
    fix_terminals();
    // TESTING this, maybe we have to compute nullable set; FIRST set computation might be different
    compute_nullable();
    // compute FIRST sets
    compute_FIRST();

    // for each non-terminal in the order in which they appear
    for (const string & nt : non_terminals) {
        cout << "FIRST(" << nt << ") = { ";
        bool firstElement = true;
        // to print the elements in grammar order, iterate over the terminals vector
        for (const string & term : terminals) {
            if (first_sets[nt].find(term) != first_sets[nt].end()) {
                if (!firstElement)
                    cout << ", ";
                cout << term;
                firstElement = false;
            }
        }
        cout << " }" << endl;
    }
}



// Task 4: FOLLOW sets
void compute_FOLLOW() {
    // ensure the classification is fixed
    fix_terminals();

    //initialize FOLLOW sets for every non-terminal
    for (const string &nt : non_terminals) {
        follow_sets[nt] = {};
    }
    // for the start symbol aka the first non terminal that shows up, add EOF marker "$"
    if (!non_terminals.empty()) {
        follow_sets[non_terminals[0]].insert("$");
    }

    bool changed;
    do {
        changed = false;
        //process each production A -> alpha
        for (const auto &rule : grammar_rules) {
            string A = rule.lhs;
            const vector<string>& alpha = rule.rhs;
            //iterate over the production's RHS
            for (size_t i = 0; i < alpha.size(); i++) {
                string B = alpha[i];
                //only process if B is a non-terminal
                if (find(non_terminals.begin(), non_terminals.end(), B) != non_terminals.end()) {
                    bool betaNullable = true;

                    //go through all the values int the alpha variable, and then if it's a terminal, then add it to the follow sets
                    for (size_t j = i + 1; j < alpha.size(); j++) {
                        string X = alpha[j];
                        //if X is a terminal, add it to FOLLOW(B) and stop.
                        if (find(terminals.begin(), terminals.end(), X) != terminals.end()) {
                            if (follow_sets[B].find(X) == follow_sets[B].end()) {
                                follow_sets[B].insert(X);
                                changed = true;

                            }
                            betaNullable = false;
                            break;
                        }
                        // if X is a non-terminal
                        if (find(non_terminals.begin(), non_terminals.end(), X) != non_terminals.end()) {
                            //add FIRST(X) (which contains only terminals) to FOLLOW(B).
                            for (const string &t : first_sets[X]) {
                                if (follow_sets[B].find(t) == follow_sets[B].end()) {
                                    follow_sets[B].insert(t);
                                    changed = true;

                                }

                            }
                            //if X is not nullable, then beta is not nullable.
                            if (nullable_set.find(X) == nullable_set.end()) {
                                betaNullable = false;
                                break;

                            }
                        }
                    }
                    // jif beta is empty or every symbol in beta is nullable, add FOLLOW(A) to FOLLOW(B)
                    if (betaNullable) {
                        //go through the folow sets and then find the union between the two
                        for (const string &t : follow_sets[A]) {
                            if (follow_sets[B].find(t) == follow_sets[B].end()) {

                                follow_sets[B].insert(t);
                                changed = true;

                            }
                        }
                    }
                }
            }
        }
    } while (changed);
}


void Task4() {
    //do the same thing as task 3, compute all the other functions from the other helper functions
    fix_terminals();
    compute_nullable();
    compute_FIRST();
    compute_FOLLOW();

    //for each non-terminal in the order of appearance, output its FOLLOW set
    for (const string &nt : non_terminals) {
        cout << "FOLLOW(" << nt << ") = { ";
        bool firstPrinted = false;
        //if $ is in the FOLLOW set then p output it first
        if (follow_sets[nt].find("$") != follow_sets[nt].end()) {
            cout << "$";
            firstPrinted = true;
        }
        //then, iterate over the terminals in the order they appear in the grammar 
        for (const string &term : terminals) {
            //skip if term is "$" or not in FOLLOW set
            if (term == "$") continue;
            if (follow_sets[nt].find(term) != follow_sets[nt].end()) {
                if (firstPrinted)
                    cout << ", " << term;
                else {
                    cout << term;
                    firstPrinted = true;
                }
            }
        }
        cout << " }" << endl;
    }
}

// Task 5: left factoring

//first compute longest common prefix between two sequences of tokens, it should take 2 vectors and return 1
vector<string> longestCommonPrefix(const vector<string>& a, const vector<string>& b) {
    vector<string> prefix;
    size_t n = min(a.size(), b.size());
    //go thorugh and if they are the same, then return it
    for (size_t i = 0; i < n; i++) {
        if (a[i] == b[i])
            prefix.push_back(a[i]);
        else
            break;
    }
    return prefix;
}

//check if vector v starts with prefix, builds off of the longest common prefix
bool startsWith(const vector<string>& v, const vector<string>& prefix) {
    if (prefix.size() > v.size()) return false;
    for (size_t i = 0; i < prefix.size(); i++) {
        if (v[i] != prefix[i])
            return false;
    }
    return true;

}

//lexicographical compare for vectors of strings
bool lexCompare(const vector<string>& a, const vector<string>& b) {
    //use the function that i found does the comparison for me
    return lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());

}

//concatenate two vectors of strings
vector<string> concat(const vector<string>& a, const vector<string>& b) {
    vector<string> result = a;
    //input the concatenated result
    result.insert(result.end(), b.begin(), b.end());
    return result;

}

//generate a new non-terminal name based on a given base
//uses a static map to ensure that if multiple new names are needed for the same base,
// the first is base1 then base2 and so on
string getNewNonTerminal(const string &base) {
    static unordered_map<string, int> new_nt_counter;
    int count = new_nt_counter[base] + 1;
    new_nt_counter[base] = count;
    return base + to_string(count);
}


//this function transforms the grammar_rules into a left-factored grammar
// sorts the resulting grammar lexicographically using the funciton above and prints the rules in the required format
void Task5() {
    //we repeatedly perform left factoring until no more factoring is possible
    bool didFactor;
    do {
        didFactor = false;
        vector<GrammarRule> newRules;
        //group productions by their lhs
        unordered_map<string, vector<GrammarRule>> groups;
        for (auto &rule : grammar_rules) {
            groups[rule.lhs].push_back(rule);
        }
        //go through and process each non-terminal group
        for (auto &pair : groups) {
            string A = pair.first;
            vector<GrammarRule>& rules = pair.second;
            //if only one production exists, no left factoring is possible
            if (rules.size() < 2) {
                for (auto &rule : rules)
                    newRules.push_back(rule);
                continue;
            }
            //find the longest common prefix among any two productions for a given value A
            vector<string> bestPrefix;
            bool foundCommon = false;
            for (size_t i = 0; i < rules.size(); i++) {
                for (size_t j = i + 1; j < rules.size(); j++) {
                    vector<string> cp = longestCommonPrefix(rules[i].rhs, rules[j].rhs);
                    if (!cp.empty()) {
                        if (!foundCommon ||
                            cp.size() > bestPrefix.size() ||
                            (cp.size() == bestPrefix.size() && lexCompare(cp, bestPrefix))) {
                            bestPrefix = cp;
                            foundCommon = true;
                        }
                    }
                }
            }
            // if no common prefix was found for A, output its productions asit was
            if (!foundCommon) {
                for (auto &rule : rules)
                    newRules.push_back(rule);
            } else {
                didFactor = true;
                //new non terminal for the left factored portion.
                string newNT = getNewNonTerminal(A);
                //partition rules into those that share the common prefix and those that do not.
                vector<GrammarRule> factorGroup;
                vector<GrammarRule> remainGroup;
                for (auto &rule : rules) {
                    if (startsWith(rule.rhs, bestPrefix))
                        factorGroup.push_back(rule);
                    else
                        remainGroup.push_back(rule);
                }
                //for the factored productions, create new productions for newNT
                vector<GrammarRule> newNTRules;
                for (auto &rule : factorGroup) {
                    vector<string> remainder(rule.rhs.begin() + bestPrefix.size(), rule.rhs.end());
                    //if remainder is empty, it represents epsilon or empty string
                    newNTRules.push_back({newNT, remainder});
                }
                //eplace the factored productions with a single production:
                // A -> bestPrefix newNT will be the grammar rule
                vector<string> newRhs = bestPrefix;
                newRhs.push_back(newNT);
                newRules.push_back({A, newRhs});
                //add back the productions that didn't factor
                for (auto &rule : remainGroup)
                    newRules.push_back(rule);
                //add the new non-terminal's productions
                for (auto &rule : newNTRules)
                    newRules.push_back(rule);
            }
        }
        grammar_rules = newRules;
    } while (didFactor);

    //sort the resulting grammar lexicographically
    sort(grammar_rules.begin(), grammar_rules.end(), [](const GrammarRule &a, const GrammarRule &b) {
        //construct sequences for lexicographical comparison
        vector<string> seqA, seqB;
        seqA.push_back(a.lhs);
        seqB.push_back(b.lhs);
        for (auto &tok : a.rhs) seqA.push_back(tok);
        for (auto &tok : b.rhs) seqB.push_back(tok);
        return lexCompare(seqA, seqB);
    });

    //print each rule in the required format:
    // LHS -> <rhs tokens separated by space> #
    for (auto &rule : grammar_rules) {
        cout << rule.lhs << " ->";
        //if the right hand side is not empty, then output the rules right hand side sequentially
        if (!rule.rhs.empty()) {
            for (auto &sym : rule.rhs)
                cout << " " << sym;
        }
        cout << " #" << endl;
    }
}


// Task 6: eliminate left recursion

string prodToString(const vector<string>& prod) {
    string s;
    for (const auto &token : prod) {
        s += token + " ";
    }
    return s;
}

void Task6() {
    // first we group productions by non terminals
    unordered_map<string, vector<vector<string>>> rules;
    unordered_set<string> originalNT;
    for (const auto &rule : grammar_rules) {
        originalNT.insert(rule.lhs);
        rules[rule.lhs].push_back(rule.rhs);
    }
    
    //get sorted list of original non-terminals which is dictionary order
    vector<string> nonterminalsSorted(originalNT.begin(), originalNT.end());
    sort(nonterminalsSorted.begin(), nonterminalsSorted.end());
    
    //for i = 1 to n over original non terminals do
    for (size_t i = 0; i < nonterminalsSorted.size(); i++) {
        string Ai = nonterminalsSorted[i];

        //for every non-terminal Aj that comes before Ai
        for (size_t j = 0; j < i; j++) {
            string Aj = nonterminalsSorted[j];

            vector<vector<string>> newProductions;
            unordered_set<string> unique;

            // For each production of Ai:
            for (auto &prod : rules[Ai]) {
                //if production is of the form Ai -> Aj γ then substitute.
                if (!prod.empty() && prod[0] == Aj) {
                    vector<string> gamma(prod.begin() + 1, prod.end());

                    //for every production for Aj (say, Aj -> δ), add Ai -> δ γ.
                    for (auto &delta : rules[Aj]) {
                        vector<string> newProd = delta; // δ
                        newProd.insert(newProd.end(), gamma.begin(), gamma.end());
                        string key = prodToString(newProd);

                        if (unique.find(key) == unique.end()) {
                            unique.insert(key);
                            newProductions.push_back(newProd);
                        }
                    }
                } else {
                    string key = prodToString(prod);
                    if (unique.find(key) == unique.end()) {
                        unique.insert(key);
                        newProductions.push_back(prod);
                    }
                }
            }
            rules[Ai] = newProductions;

        }
        
        //eliminate immediate left recursion from productions of Ai 
        vector<vector<string>> leftRecursive;
        vector<vector<string>> nonLeftRecursive;
        for (auto &prod : rules[Ai]) {
            if (!prod.empty() && prod[0] == Ai)
                leftRecursive.push_back(vector<string>(prod.begin() + 1, prod.end()));
            else
                nonLeftRecursive.push_back(prod);
        }
        
        if (!leftRecursive.empty()) {
            //introduce a new non-terminal name following the convention: Ai1
            string Ai1 = Ai + "1";
            vector<vector<string>> newAiProductions;
            //for each non-left-recursive production b, replace it with: Ai -> b Ai1
            {
                unordered_set<string> unique;
                for (auto &beta : nonLeftRecursive) {
                    vector<string> newBeta = beta;
                    newBeta.push_back(Ai1);
                    string key = prodToString(newBeta);

                    if (unique.find(key) == unique.end()) {
                        unique.insert(key);
                        newAiProductions.push_back(newBeta);

                    }
                }
            }
            rules[Ai] = newAiProductions;
            
            //for each left-recursive production of the form Ai -> Ai a, lets create a production for Ai1: Ai1 -> a Ai1 
            vector<vector<string>> Ai1Productions;
            {
                unordered_set<string> unique;
                for (auto &alpha : leftRecursive) {
                    vector<string> newAlpha = alpha;
                    newAlpha.push_back(Ai1);
                    string key = prodToString(newAlpha);
                    if (unique.find(key) == unique.end()) {
                        unique.insert(key);
                        Ai1Productions.push_back(newAlpha);
                    }
                }
                //add the epsilon or empty production represented as an empty vector
                string epsKey = prodToString(vector<string>());
                if (unique.find(epsKey) == unique.end()) {
                    unique.insert(epsKey);
                    Ai1Productions.push_back(vector<string>());
                }
            }
            //add the new non-terminal's productions to the rules
            rules[Ai1] = Ai1Productions;

        }
    }
    
    //collect all rules into a vector of GrammarRule
    vector<GrammarRule> finalRules;
    
    for (auto &entry : rules) {
        const string &lhs = entry.first;

        for (auto &prod : entry.second) {
            finalRules.push_back({lhs, prod});
        }
    }
    
    //sort the final rules lexicographically so call the function we defined above
    sort(finalRules.begin(), finalRules.end(), [](const GrammarRule &a, const GrammarRule &b) {
        if (a.lhs != b.lhs)
            return a.lhs < b.lhs;

        return lexicographical_compare(a.rhs.begin(), a.rhs.end(), b.rhs.begin(), b.rhs.end());

    });
    
    //print the final grammar
    //each rule is printed as: LHS -> <rhs symbols separated by a space> #
    //if the RHS is empty or epsilon, print nothing between -> and #
    for (auto &rule : finalRules) {
        cout << rule.lhs << " ->";
        if (!rule.rhs.empty()) {
            for (auto &sym : rule.rhs)
                cout << " " << sym;
        }
        cout << " #" << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Error: missing argument\n";
        return 1;
    }

    
    /*
       Note that by convention argv[0] is the name of your executable,
       and the first argument to your program is stored in argv[1]
     */

    int task = atoi(argv[1]);


    ReadGrammar();  // Reads the input grammar from standard input
                    // and represent it internally in data structures
                    // as described in project 2 presentation file
    

    switch (task) {
        case 1: Task1(); 
            break;
        case 2: Task2(); 
            break;
        case 3: Task3(); 
            break;
        case 4: Task4(); 
            break;
        case 5: Task5(); 
            break;
        case 6: Task6(); 
            break;
        default:
            cout << "Error: unrecognized task number " << task << "\n";
            break;
    }

    return 0;
}