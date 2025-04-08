# DAA!
## Bellman Ford
```cpp
#include <iostream>
#include <vector>
using namespace std;

struct Edge {
    int src;
    int dest;
    int wt;
};

vector<int> bellman(vector<Edge>& edges, int V, int E, int S) {
    vector<int> dist(V, 1e8);
    dist[S] = 0;
  
    for (int i=0; i<V-1; i++) {
        for (auto edge: edges) {
            int u = edge.src;
            int v = edge.dest;
            int wt = edge.wt;
  
            if (dist[u] != 1e8 && dist[u] + wt < dist[v]) {
                dist[v] = dist[u] + wt;
            }
        }
    }
  
    for (auto edge: edges) {
        int u = edge.src;
        int v = edge.dest;
        int wt = edge.wt;
        if (dist[u] != 1e8 && dist[u] + wt < dist[v]) {
            return {-1};
        }
    }
    return dist;
}
  
int main() {
    int V, E, S;
    cout << "Enter number of vertices: ";
    cin >> V;
    cout << "Enter number of edges: ";
    cin >> E;
    vector<Edge> edges(E);
    cout << "Enter edge details (source destination weight):" << endl;
    for (int i = 0; i < E; i++) {
        cin >> edges[i].src >> edges[i].dest >> edges[i].wt;
    }
    cout << "Enter source vertex: ";
    cin >> S;
    vector<int> distances = bellman(edges, V, E, S);
    if (distances.size() == 1 && distances[0] == -1) {
        cout << "Negative cycle detected in the graph." << endl;
    } else {
        cout << "Shortest distances from source vertex " << S << ":" << endl;
        for (int i = 0; i < V; i++) {
            if (distances[i] == 1e8) {
                cout << "Vertex " << i << ": INFINITY" << endl;
            } else {
                cout << "Vertex " << i << ": " << distances[i] << endl;
            }
        }
    }
    return 0;
}
```

## Ford-Fulkerson; Edmond-Karp
```cpp
#include <iostream>
#include <vector>
#include <stack>
#include <queue>
using namespace std;

int dfs(vector<vector<int>>& graph, int src, int sink, vector<int>& parent) {
    int n = graph.size();
    vector<int> visited(n, 0);
    stack<int> s;

    visited[src] = 1;
    parent[src] = -1;
    s.push(src);

    while (!s.empty()) {
        int u = s.top();
        s.pop();

        for (int v=0; v<n; v++) {
            if (!visited[v] && graph[u][v]>0) {
                visited[v] = 1;
                parent[v] = u;
                s.push(v);
            }
        }
    }

    return visited[sink];
}

// BFS function for Edmonds-Karp algorithm
int bfs(vector<vector<int>>& graph, int src, int sink, vector<int>& parent) {
    int n = graph.size();
    fill(parent.begin(), parent.end(), -1);
    vector<int> visited(n, 0);
    queue<int> q;
    
    visited[src] = 1;
    parent[src] = -1;
    q.push(src);
    
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        
        for (int v = 0; v < n; v++) {
            if (!visited[v] && graph[u][v] > 0) {
                visited[v] = 1;
                parent[v] = u;
                q.push(v);
            }
        }
    }
    
    return visited[sink];
}

// Helper function to print a single path
void printPath(vector<int>& parent, int src, int sink) {
    vector<int> path;
    for (int v = sink; v != src; v = parent[v]) {
        path.push_back(v);
    }
    path.push_back(src);
    
    cout << "Path: ";
    for (int i = path.size() - 1; i >= 0; i--) {
        cout << path[i];
        if (i > 0) cout << " -> ";
    }
}

// Modified FF function to track and print paths
int FF(vector<vector<int>>& graph, int src, int sink, vector<int>& parent) {
    int n = graph.size();
    int maxFlow = 0;
    int pathCount = 0;

    while (dfs(graph, src, sink, parent)) {
        int pathFlow = INT_MAX;
        for (int v=sink; v!=src; v = parent[v]) {
            int u = parent[v];
            pathFlow = min(pathFlow, graph[u][v]);
        }

        // Print the current augmenting path and its flow
        cout << "Augmenting Path #" << ++pathCount << " (Flow: " << pathFlow << "): ";
        printPath(parent, src, sink);
        cout << endl;

        for (int v=sink; v!=src; v=parent[v]) {
            int u = parent[v];
            graph[u][v] -= pathFlow;
            graph[v][u] += pathFlow;
        }

        maxFlow += pathFlow;
    }
    return maxFlow;
}

// Edmonds-Karp algorithm (Ford-Fulkerson with BFS)
int EK(vector<vector<int>>& graph, int src, int sink, vector<int>& parent) {
    int n = graph.size();
    int maxFlow = 0;
    int pathCount = 0;

    while (bfs(graph, src, sink, parent)) {
        int pathFlow = INT_MAX;
        for (int v = sink; v != src; v = parent[v]) {
            int u = parent[v];
            pathFlow = min(pathFlow, graph[u][v]);
        }

        // Print the current augmenting path and its flow
        cout << "Augmenting Path #" << ++pathCount << " (Flow: " << pathFlow << "): ";
        printPath(parent, src, sink);
        cout << endl;

        for (int v = sink; v != src; v = parent[v]) {
            int u = parent[v];
            graph[u][v] -= pathFlow;
            graph[v][u] += pathFlow;
        }

        maxFlow += pathFlow;
    }
    return maxFlow;
}

// Function to print all flow paths and residual graph
void printFlowPaths(vector<vector<int>>& residualGraph, int n) {
    cout << "\nResidual Graph:" << endl;
    cout << "  ";
    for (int i = 0; i < n; i++) {
        cout << i << " ";
    }
    cout << endl;
    
    for (int i = 0; i < n; i++) {
        cout << i << " ";
        for (int j = 0; j < n; j++) {
            cout << residualGraph[i][j] << " ";
        }
        cout << endl;
    }
}

int main() {
    int N, E;
    cin >> N >> E;

    vector<vector<int>> graph(N, vector<int>(N, 0));
    vector<int> parent(N, -1);

    for (int i = 0; i < E; i++) {
        int u, v, c;
        cin >> u >> v >> c;
        graph[u][v] = c;
    }

    int src = 0, sink = N - 1;
    
    cout << "Choose algorithm:\n";
    cout << "1. Ford-Fulkerson with DFS\n";
    cout << "2. Edmonds-Karp with BFS\n";
    int choice;
    cin >> choice;
    
    int maxFlow;
    if (choice == 2) {
        cout << "\nRunning Edmonds-Karp Algorithm...\n";
        maxFlow = EK(graph, src, sink, parent);
    } else {
        cout << "\nRunning Ford-Fulkerson Algorithm with DFS...\n";
        maxFlow = FF(graph, src, sink, parent);
    }
    
    cout << "\nMaximum Flow: " << maxFlow << endl;
    
    // Print the residual graph after finding max flow
    printFlowPaths(graph, N);

    return 0;
}
```

## KMP
```cpp
#include <iostream>
#include <vector>
using namespace std;

vector<int> computeLPS(string pat) {
    int m = pat.size();
    vector<int> lps(m, 0);
    int prevLPS = 0;

    for (int i=1; i<m;) {
        if (pat[i] == pat[prevLPS]) {
            lps[i++] = ++prevLPS;
        } else if (prevLPS) {
            prevLPS = lps[prevLPS-1];
        } else {
            lps[i++] = 0;
        }
    }
    return lps;
}

void KMP(string txt, string pat) {
    int n = txt.size(), m = pat.size();
    bool found = false;
    vector<int> lps = computeLPS(pat);

    for (int i=0, j=0; i<n;) {
        if (txt[i] == pat[j]) {
            i++;
            j++;
        } else if (i<n && txt[i] != pat[j]) {
            if (j) {
                j = lps[j-1];
            } else {
                i++;
            }
        }

        if (j == m) {
            cout << "Found pattern at index " << i-j << endl;
            found = true;
            j = lps[j-1];
        }
    }

    if (!found) {
        cout<<"Not found"<<endl;
    }
}

int main() {
    string txt, pat;
    cin>>txt;
    cin>>pat;
    KMP(txt, pat);
    return 0;
}
```

## MCM
```cpp
#include <iostream>
#include <vector>
using namespace std;

int lookup(vector<vector<int>>& m, vector<vector<int>>& s, vector<int>& p, int i, int j) {
    if (i==j) return 0;
    if (m[i][j] != INT_MAX) return m[i][j];

    for (int k=i; k<j; k++) {
        int cost = lookup(m, s, p, i, k) + lookup(m, s, p, k+1, j) + p[i-1]*p[k]*p[j];
        if (cost < m[i][j]) {
            m[i][j] = cost;
            s[i][j] = k;
        }
    }

    return m[i][j];
}

void print_prarantheses(vector<vector<int>>& s, int i, int j) {
    if (i==j) {
        cout<<'A'<<i;
    }
    else {
        cout<<'(';
        print_prarantheses(s, i, s[i][j]);
        print_prarantheses(s, s[i][j]+1, j);
        cout<<')';
    }
}

int main() {
    vector<int> p;
    while (true) {
        int temp;
        cin>>temp;
        if (temp==-1) break;
        p.push_back(temp);
    }
    int n = p.size()-1;

    vector<vector<int>> m(n+1, vector<int>(n+1, INT_MAX));
    vector<vector<int>> s(n+1, vector<int>(n+1, 0));

    int minCost = lookup(m, s, p, 1, n);
    cout<<minCost<<endl;
    print_prarantheses(s, 1, n);

    return 0;
}
```

## N-Queens
```cpp
#include <iostream>
#include <vector>
#include <unordered_set>
using namespace std;

void backtrack(
    int r, int n,
    unordered_set<int>& col,
    unordered_set<int>& posDiag,
    unordered_set<int>& negDiag,
    vector<string>& board,
    vector<vector<string>>& res
) {
    if (r == n) {
        res.push_back(board);
        return;
    }

    for (int c=0; c<n; c++) {
        if (col.count(c) || posDiag.count(r+c) || negDiag.count(r-c)) continue;

        col.insert(c);
        posDiag.insert(r+c);
        negDiag.insert(r-c);
        board[r][c] = 'Q';

        backtrack(r+1, n, col, posDiag, negDiag, board, res);

        col.erase(c);
        posDiag.erase(r+c);
        negDiag.erase(r-c);
        board[r][c] = '.';

    }
}


int main() {
    int n = 8;
    unordered_set<int> col, posDiag, negDiag;
    vector<string> board(n, string(n, '.'));
    vector<vector<string>> res;

    backtrack(0, n, col, posDiag, negDiag, board, res);

    cout << "Found " << res.size() << " solutions for " << n << "-Queens problem:" << endl << endl;

    for (int i = 0; i < res.size(); i++) {
        cout << "Solution " << (i+1) << ":" << endl;
        for (const auto& row: res[i]) {
            cout << row << endl;
        }
        cout << endl;
    }
    
    if (res.empty()) {
        cout << "No solutions found!" << endl;
    }
    return 0;
}

```

## Naive String Matching
```cpp
#include <iostream>
#include <vector>
using namespace std;

void NSM(string txt, string pat) {
    bool found = false;
    int n = txt.size(), m = pat.size();
    
    for (int i=0; i<n-m+1; i++) {
        int j;
        for (j=0; j<m; j++) {
            if (txt[i+j] != pat[j]) {
                break;
            }
        }
        if (j==m) {
            cout<<"Pattern found at "<<i<<endl;
            found = true;
        }
    }
    
    if (!found) {
        cout<<"Pattern not found"<<endl;
    }
}

int main() {
    string txt, pat;
    cin>>txt;
    cin>>pat;
    NSM(txt, pat);
    return 0;
}
```

## Max Subarray
```cpp
#include<iostream>
using namespace std;
void maxSubarry(int *arr,int n){
    int max_sum=arr[0];
    int curr_sum=arr[0];
    int start=0, end=0, s=0;
    for(int i=0;i<n;i++){
        curr_sum+=arr[i];
        if(curr_sum>max_sum){
            max_sum=curr_sum;
            start=s;
            end=i;
        }
        if(curr_sum<=0){
            curr_sum=0;
            s=i+1;
        }
    }
    cout<< "Maximum sum: "<<max_sum<<endl;
    cout<< "Start index: "<<start<<", "<<"End Index: "<<end<<endl;
}
int main(){
    int n;
    cin>>n;
    int* arr=new int[n];
    for(int i=0;i<n;i++){
        cin>>arr[i];
    }
    maxSubarry(arr,n);
    return 0;
}
```
## Karatsuba
```cpp
#include<iostream>
#include<cmath>
using namespace std;
long long karatsuba(long long a, long long b) {
    if(a<10 || b <10){
        return a*b;
    }
    int n=max((int)log10(a)+1,(int)log10(b)+1);
    int m=n/2;
    long long mul=pow(10,m);
    long long ar=a%mul;
    long long br=b%mul;
    long long al=a/mul;
    long long bl=b/mul;
    long long z0=karatsuba(al,bl);
    long long z1=karatsuba(al+ar,bl+br);
    long long z2=karatsuba(ar,br);
    return (pow(10,n)*z2+pow(10,m)*(z1-z2-z0)+z0);
}
int main(){
    int x,y;
    cin>>x>>y;
    cout<<karatsuba(x,y);
    return 0;
}
```

## LCS
```cpp
#include<iostream>
#include<vector>
#include<string>
#include<set>
using namespace std;
set<string> LCS(string a, string b){
    int m=a.size(), n=b.size();
    vector<vector<int>> dp(m+1,vector<int>(n+1,0));
    vector<vector<set<string>>> lcsres(m+1,vector<set<string>>(n+1));
    for(int i=0;i<=m;i++){
        for(int j=0;j<=n;j++){
            if(i==0 || j==0){
                lcsres[i][j].insert("");
            } else if (a[i-1]==b[j-1]) {
                for (const string& s : lcsres[i - 1][j - 1])
                    lcsres[i][j].insert(s + a[i - 1]);
                dp[i][j] = dp[i - 1][j - 1] + 1;
            } else {
                if (dp[i - 1][j] >= dp[i][j - 1]) {
                    dp[i][j] = dp[i - 1][j];
                    lcsres[i][j].insert(lcsres[i-1][j].begin(),lcsres[i-1][j].end());
                } if (dp[i - 1][j] <= dp[i][j - 1]) {
                    dp[i][j] = dp[i][j - 1];
                    lcsres[i][j].insert(lcsres[i][j - 1].begin(), lcsres[i][j - 1].end());
                }
            }
        }
    }
    return lcsres[m][n];
}
int main() {
    string x, y;
    getline(cin, x);
    getline(cin, y);
    auto result = LCS(x, y);
    cout << "The length of the Longest Common Subsequence is: " << (*result.begin()).size() << endl;
    for (const string& s : result)
        cout << s << endl;
    return 0;
}
```

## KMP
```cpp
#include<iostream>
#include<vector>
#include<string>
using namespace std;
vector<int> computeLPS(string pat){
    int m=pat.size(),len=0;
    vector<int> lps(m,0);
    for(int i=1;i<m;){
        if(pat[i]==pat[len]){
            lps[i++]=++len;
        } else if(len){
            len=lps[len-1];
        } else {
            lps[i++]=0;
        }
    }
    return lps;
}
void KMPSearch(string txt, string pat){
    int n=txt.size(),m=pat.size();
    bool found = false;
    vector<int> LPS = computeLPS(pat);
    for(int i=0,j=0;i<n;){
        if(txt[i]==pat[j]){
            i++;
            j++;
        }
        if(j==m){
            cout << "Found pattern at index " << i-j << endl;
            found=true;
            j=LPS[j-1];
        } else if (i<n && txt[i]!=pat[j]){
            if (j) j=LPS[j-1];
            else i++;
        }
    }
    if (!found) cout << "Pattern not found" << endl;
}
int main(){
    string a,b;
    getline(cin,a);
    getline(cin,b);
    KMPSearch(a, b);
    return 0;
}
```

## Rabin Karp
```cpp
#include<iostream>
#include<string>
using namespace std;
void RabinKarp(string txt, string pat){
    int n=txt.size(),m=pat.size();
    int h=1; int p=0,t=0;
    for(int i=0;i<m-1;i++) h=(h*10)%13;
    for(int i=0;i<m;i++) {
        p=(10*p+pat[i])%13;
        t=(10*t+txt[i])%13;
    }
    for(int i=0;i<n-m;i++){
        if(p==t){
            if(txt.substr(i,m)==pat) cout<<i<<" ";
        }
        if(i<n-m){
            t=(10*(t-txt[i]*h)+txt[i+m])%13;
            if(t<0) t+=13;
        }
    }
}
int main() {
    string txt, pat;
    getline(cin, txt);
    getline(cin, pat);
    RabinKarp(txt, pat);
    return 0;
}
```

## MCM
```cpp
#include<iostream>
#include<vector>
using namespace std;
void printmat(vector<vector<int>> &s,int i, int j){
    if(i==j){
        cout<<'A'<<i;
        return;
    }
    cout<<'(';
    printmat(s,i,s[i][j]);
    printmat(s,s[i][j]+1,j);
    cout<<')';
}
int main(){
    vector<int> p;
    int x;
    while(cin>>x && x!=-1) p.push_back(x);
    int n=p.size()-1;
    vector<vector<int>> m(n+1,vector<int>(n+1,0));
    vector<vector<int>> s(n+1,vector<int>(n+1,0));
    for(int len=2;len<=n;len++){
        for(int i=1;i<=n-len+1;i++){
            int j=i+len-1;
            m[i][j] = 1e9;
            for(int k=i;k<j;k++){
                int c = m[i][k]+m[k+1][j]+p[i-1]*p[k]*p[j];
                if(c< m[i][j]){
                    m[i][j]=c;
                    s[i][j]=k;
                }
            }
        }
    }
    cout << m[1][n] << endl;
    printmat(s, 1, n);
    return 0;
}
```
## Flloyd Warshall
```cpp
#include<iostream>
#include<vector>
using namespace std;
int main(){
    int n;
    cin>>n;
    vector<vector<int>> dist(n,vector<int>(n,0));
    vector<vector<int>> next(n,vector<int>(n,-1));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            cin >> dist[i][j];
            if(dist[i][j]!=999 && i!=j) next[i][j]=j;
        }
    }
    for(int k=0;k<n;k++){
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                if(dist[i][j]>dist[i][k]+dist[k][j]){
                    dist[i][j]=dist[i][k]+dist[k][j];
                    next[i][j]=next[i][k];
                }
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            cout << dist[i][j] << " ";
        cout << endl;
    }
    cout << "Shortest paths:\n";
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (i != j) {
                cout << i << " to " << j << ": " << i;
                for (int k = next[i][j]; k!=-1 && k != j; k = next[k][j])
                    cout << " -> " << k;
                cout << " -> " << j << endl;
            }
        }
    }   
    return 0;
}
```

## Jarvis
```cpp
#include<iostream>
#include<vector>
using namespace std;
struct Point
{
    int x,y;
};
int orient(Point a, Point b, Point c){
    int val=(b.y-a.y)*(c.x-b.x)-(b.x-a.x)*(c.y-b.y);
    return val==0?0:(val<0?1:2);
}
int convexHull(vector<Point>& points){
    int n=points.size();
    if (n<3) return n;
    int start=0;
    for(int i=1;i<n;i++){
        if(points[i].x<points[start].x){
            start=i;
        }
    }
    vector<Point> hull;
    int p=start;
    do {
        hull.push_back(points[p]);
        int q=(p+1)%n;
        for(int i=0;i<n;i++){
            if(orient(points[p],points[i],points[q])==2) q=i;
        }
        p=q;
    }while(p!=start);
    return hull.size();
}
int main(){
    int n;
    cin >> n;
    vector<Point> points(n);
    for (int i = 0; i < n; i++)
        cin >> points[i].x >> points[i].y;
    cout << convexHull(points);
    return 0;
}
```

## N Queens
```cpp
#include<iostream>
#include<vector>
#include<unordered_set>
using namespace std;
void solve(int row, int n, 
    unordered_set<int> &cols, 
    unordered_set<int> &pdiag, 
    unordered_set<int> &ndiag, 
    vector<string> &board,
    vector<vector<string>> &result 
) {
    if(row==n){
        result.push_back(board);
        return;
    }
    for(int col=0;col<n;col++){
        if(cols.count(col)||pdiag.count(row+col)||ndiag.count(row-col)){
            continue;
        }
        cols.insert(col);
        pdiag.insert(row+col);
        ndiag.insert(row-col);
        board[row][col]='Q';
        solve(row+1,n,cols,pdiag,ndiag,board,result);
        board[row][col]='.';
        cols.erase(col);
        pdiag.erase(row+col);
        ndiag.erase(row-col);
    }
}
int main(){
    int n=8;
    vector<vector<string>> result;
    vector<string> board(n,string(n,'.'));
    unordered_set<int> cols,pdiag,ndiag;
    solve(0,n,cols,pdiag,ndiag,board, result);
    cout<<"Solution Found: "<<result.size()<<endl;
    for(auto& sol: result){
        for(auto& row: sol){
            cout<<row<<endl;
        }
        cout<<endl;
    }
    return 0;
}
```

## Job Selection
```cpp
#include<iostream>
#include<vector>
#include<algorithm>
using namespace std;
struct Job {
    int profit,deadline,resources;
};
void explore(vector<Job>&jobs,int M,int i, int profit, int resources, int &maxProfit, string indent=""){
    cout << indent << "Job " << i << ": " 
         << (i == 0 ? "Start" : (resources + jobs[i-1].resources <= M ? "Included" : "Excluded")) 
         << " | Profit: " << profit 
         << " | Resources Used: " << resources << endl;
    if(i==jobs.size()){
        maxProfit=max(maxProfit,profit);
        return;
    }
    if(jobs[i].resources+resources<=M){
        explore(jobs,M,i+1,profit+jobs[i].profit, resources+jobs[i].resources,maxProfit, indent+" ");
    }
    explore(jobs,M,i+1,profit,resources,maxProfit,indent+" ");
}
int main(){
    int n,M;
    cin>>n>>M;
    vector<Job> jobs(n);
    for(int i=0;i<n;i++){
        cin>>jobs[i].profit>>jobs[i].deadline>>jobs[i].resources;
    }
    sort(jobs.begin(), jobs.end(),[](const Job&a, const Job&b){
        return a.profit>b.profit;
    });
    int maxProfit=0;
    explore(jobs,M,0,0,0,maxProfit);
    cout<<maxProfit;
    return 0;
}
```
## Randomized Quicksort
```cpp
#include<iostream>
#include<vector>
#include<ctime>
using namespace std;
void qsort(vector<int> &a,int l,int h){
    if(l>=h) return;
    int r=l+rand()%(h-l+1);
    swap(a[r],a[h]);
    int p=a[h], i=l;
    for(int j=l;j<h;j++){
        if(a[j]<=p) swap(a[i++],a[j]);
    }
    swap(a[i],a[h]);
    qsort(a,l,i-1);
    qsort(a,i+1,h);
}
int main(){
    srand(time(0));
    int n; cin>>n;
    vector<int> a(n);
    for(int i=0;i<n;i++) cin>>a[i];
    qsort(a,0,n-1);
    for(int x:a) cout<<x<<" ";
    return 0;
}
```

## Cheat sort
```cpp
#include<iostream>
#include<vector>
#include<algorithm>
using namespace std;
int main(){
    int n; cin>>n;
    vector<int> nums(n,0);
    for(int i=0;i<n;i++) cin>>nums[i];
    sort(nums.begin(),nums.end());
    for(int x:nums) cout<<x<<" ";
    return 0;
}
```

## AI Powered Karp
```cpp
#include <iostream>
#include <string>
using namespace std;

int main() {
    string text = "The quick brown fox jumps over the lazy dog";
    string pattern = "quick";

    size_t pos = text.find(pattern);
    if (pos != string::npos) {
        cout << "Pattern found at index: " << pos << endl;
    } else {
        cout << "Pattern not found" << endl;
    }

    return 0;
}
```
