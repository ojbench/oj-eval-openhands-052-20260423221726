#include<iostream>
#include<algorithm>
#include<cstring>
#include<string>
#include<vector>

class term {
 public:
    int a, b, c, d;

    term() {}
    term(int _a, int _b, int _c, int _d) {
        a = _a; b = _b; c = _c; d = _d;
    }
    bool operator == (const term &obj) const {
        return b == obj.b && c == obj.c && d == obj.d;
    }
    bool operator != (const term &obj) const {
        return b != obj.b || c != obj.c || d != obj.d;
    }
    bool operator < (const term &obj) const {
        if (b != obj.b) return b > obj.b;
        if (c != obj.c) return c > obj.c;
        return d > obj.d;
    }
};

class poly {
 public:
    int n;
    term *t;

    poly() {n = 0; t = NULL;}
    poly(int _n) {
        n = _n;
        t = new term[n];
    }
    poly(const poly &p) {
        n = p.n;
        t = new term[n];
        for (int i = 0; i < n; ++i) {
            t[i] = p.t[i];
        }
    }
    void simplify() {
        if (n == 0) return;
        std::sort(t, t + n);
        int cnt = 0;
        for (int i = 0; i < n; ++i) {
            if (i == 0 || t[i] != t[cnt - 1]) {
                t[cnt++] = t[i];
            } else {
                t[cnt - 1].a += t[i].a;
            }
        }
        // Remove zero terms
        int new_cnt = 0;
        for (int i = 0; i < cnt; ++i) {
            if (t[i].a != 0) {
                t[new_cnt++] = t[i];
            }
        }
        n = new_cnt;
    }
    poly operator + (const poly &obj) const {
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) ans.t[i + n] = obj.t[i];
        ans.simplify();
        return ans;
    }
    poly operator - (const poly &obj) const {
        poly ans(n + obj.n);
        for (int i = 0; i < n; ++i) ans.t[i] = t[i];
        for (int i = 0; i < obj.n; ++i) {
            ans.t[i + n] = obj.t[i];
            ans.t[i + n].a *= -1;
        }
        ans.simplify();
        return ans;
    }
    poly operator * (const poly &obj) const {
        poly ans(n * obj.n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < obj.n; ++j) {
                ans.t[i * obj.n + j].a = t[i].a * obj.t[j].a;
                ans.t[i * obj.n + j].b = t[i].b + obj.t[j].b;
                ans.t[i * obj.n + j].c = t[i].c + obj.t[j].c;
                ans.t[i * obj.n + j].d = t[i].d + obj.t[j].d;
            }
        }
        ans.simplify();
        return ans;
    }
    poly& operator = (const poly &obj) {
        if (&obj == this) return *this;
        n = obj.n;
        if (t != NULL) delete []t;
        t = new term[n];
        for (int i = 0; i < n; ++i) {
            t[i] = obj.t[i];
        }
        return *this;
    }
    poly derivate() const {
        std::vector<term> result;
        for (int i = 0; i < n; ++i) {
            // Derivative of a*x^b*sin^c(x)*cos^d(x)
            // Using product rule: (uvw)' = u'vw + uv'w + uvw'
            
            // x^b derivative: b*x^(b-1)
            if (t[i].b > 0) {
                term nt;
                nt.a = t[i].a * t[i].b;
                nt.b = t[i].b - 1;
                nt.c = t[i].c;
                nt.d = t[i].d;
                result.push_back(nt);
            }
            
            // sin^c(x) derivative: c*sin^(c-1)(x)*cos(x)
            if (t[i].c > 0) {
                term nt;
                nt.a = t[i].a * t[i].c;
                nt.b = t[i].b;
                nt.c = t[i].c - 1;
                nt.d = t[i].d + 1;
                result.push_back(nt);
            }
            
            // cos^d(x) derivative: -d*sin(x)*cos^(d-1)(x)
            if (t[i].d > 0) {
                term nt;
                nt.a = -t[i].a * t[i].d;
                nt.b = t[i].b;
                nt.c = t[i].c + 1;
                nt.d = t[i].d - 1;
                result.push_back(nt);
            }
        }
        
        poly ans(result.size());
        for (size_t i = 0; i < result.size(); ++i) {
            ans.t[i] = result[i];
        }
        ans.simplify();
        return ans;
    }
    ~poly() {
        if (t != NULL) delete []t;
    }
};

class frac {
 public:
    poly p, q;

    frac() {}
    frac(int x) {
        p = poly(1);
        p.t[0] = term(x, 0, 0, 0);
        q = poly(1);
        q.t[0] = term(1, 0, 0, 0);
    }
    frac(term _p) {
        q = poly(1);
        q.t[0] = term(1, 0, 0, 0);
        p = poly(1);
        p.t[0] = _p;
    }
    frac(poly _p, poly _q) : p(_p), q(_q) {}

    frac operator + (const frac &obj) const {
        return frac(p * obj.q + q * obj.p, q * obj.q);
    }
    frac operator - (const frac &obj) const {
        return frac(p * obj.q - q * obj.p, q * obj.q);
    }
    frac operator * (const frac &obj) const {
        return frac(p * obj.p, q * obj.q);
    }
    frac operator / (const frac &obj) const {
        return frac(p * obj.q, q * obj.p);
    }
    frac derivate() const {
        // (p/q)' = (p'*q - p*q')/(q*q)
        poly p_prime = p.derivate();
        poly q_prime = q.derivate();
        return frac(p_prime * q - p * q_prime, q * q);
    }
    void output() {
        // Check if numerator is 0
        if (p.n == 0) {
            std::cout << "0" << std::endl;
            return;
        }
        
        // Check if denominator is 1
        bool denom_is_one = (q.n == 1 && q.t[0].a == 1 && q.t[0].b == 0 && q.t[0].c == 0 && q.t[0].d == 0);
        
        if (denom_is_one) {
            output_poly(p);
            std::cout << std::endl;
        } else {
            // Output with parentheses if needed
            if (p.n > 1) std::cout << "(";
            output_poly(p);
            if (p.n > 1) std::cout << ")";
            
            std::cout << "/";
            
            if (q.n > 1) std::cout << "(";
            output_poly(q);
            if (q.n > 1) std::cout << ")";
            
            std::cout << std::endl;
        }
    }
    
    void output_poly(const poly &pol) {
        for (int i = 0; i < pol.n; ++i) {
            if (i > 0) {
                if (pol.t[i].a > 0) {
                    std::cout << "+";
                }
            }
            
            int coef = pol.t[i].a;
            int xpow = pol.t[i].b;
            int sinpow = pol.t[i].c;
            int cospow = pol.t[i].d;
            
            // Check if it's a constant term
            bool is_const = (xpow == 0 && sinpow == 0 && cospow == 0);
            
            // Output coefficient
            if (is_const) {
                std::cout << coef;
            } else {
                if (coef == 1) {
                    // Don't output anything
                } else if (coef == -1) {
                    std::cout << "-";
                } else {
                    std::cout << coef;
                }
            }
            
            // Output x^b
            if (xpow > 0) {
                std::cout << "x";
                if (xpow > 1) {
                    std::cout << "^" << xpow;
                }
            }
            
            // Output sin^c(x)
            if (sinpow > 0) {
                std::cout << "sin";
                if (sinpow > 1) {
                    std::cout << "^" << sinpow;
                }
                std::cout << "x";
            }
            
            // Output cos^d(x)
            if (cospow > 0) {
                std::cout << "cos";
                if (cospow > 1) {
                    std::cout << "^" << cospow;
                }
                std::cout << "x";
            }
        }
    }
};

// Parser
int pos;
char *s;
int n;

frac parse_expr();
frac parse_term_expr();
frac parse_factor();

// Skip whitespace (if any)
void skip_space() {
    while (pos < n && s[pos] == ' ') pos++;
}

// Parse a positive number (no sign)
int parse_positive_number() {
    int num = 0;
    while (pos < n && s[pos] >= '0' && s[pos] <= '9') {
        num = num * 10 + (s[pos] - '0');
        pos++;
    }
    return num;
}

// Parse a basic term (coefficient * x^b * sin^c * cos^d)
frac parse_basic_term() {
    skip_space();
    
    int coef = 1;
    int xpow = 0, sinpow = 0, cospow = 0;
    
    // Check for coefficient
    if (s[pos] >= '0' && s[pos] <= '9') {
        coef = parse_positive_number();
    }
    
    // Parse x, sin, cos in any order
    while (pos < n) {
        if (s[pos] == 'x') {
            pos++;
            int pow = 1;
            if (pos < n && s[pos] == '^') {
                pos++;
                pow = parse_positive_number();
            }
            xpow += pow;
        } else if (pos + 2 < n && s[pos] == 's' && s[pos+1] == 'i' && s[pos+2] == 'n') {
            pos += 3;
            int pow = 1;
            if (pos < n && s[pos] == '^') {
                pos++;
                pow = parse_positive_number();
            }
            if (pos < n && s[pos] == 'x') pos++;
            sinpow += pow;
        } else if (pos + 2 < n && s[pos] == 'c' && s[pos+1] == 'o' && s[pos+2] == 's') {
            pos += 3;
            int pow = 1;
            if (pos < n && s[pos] == '^') {
                pos++;
                pow = parse_positive_number();
            }
            if (pos < n && s[pos] == 'x') pos++;
            cospow += pow;
        } else {
            break;
        }
    }
    
    term t(coef, xpow, sinpow, cospow);
    return frac(t);
}

// Parse factor (handles parentheses and basic terms)
frac parse_factor() {
    skip_space();
    
    if (s[pos] == '(') {
        pos++;
        frac result = parse_expr();
        if (pos < n && s[pos] == ')') pos++;
        return result;
    }
    
    if (s[pos] == '-') {
        pos++;
        frac result = parse_factor();
        poly zero(0);
        result.p = zero - result.p;
        return result;
    }
    
    if (s[pos] == '+') {
        pos++;
        return parse_factor();
    }
    
    return parse_basic_term();
}

// Parse multiplication and division
frac parse_term_expr() {
    frac left = parse_factor();
    
    while (pos < n) {
        skip_space();
        if (pos >= n) break;
        
        if (s[pos] == '*') {
            pos++;
            frac right = parse_factor();
            left = left * right;
        } else if (s[pos] == '/') {
            pos++;
            frac right = parse_factor();
            left = left / right;
        } else {
            break;
        }
    }
    
    return left;
}

// Parse addition and subtraction
frac parse_expr() {
    frac left = parse_term_expr();
    
    while (pos < n) {
        skip_space();
        if (pos >= n) break;
        
        if (s[pos] == '+') {
            pos++;
            frac right = parse_term_expr();
            left = left + right;
        } else if (s[pos] == '-') {
            pos++;
            frac right = parse_term_expr();
            left = left - right;
        } else {
            break;
        }
    }
    
    return left;
}

void solve(char *str, int len) {
    s = str;
    n = len;
    pos = 0;
    
    frac f = parse_expr();
    f.output();
    
    frac f_prime = f.derivate();
    f_prime.output();
}

int main() {
    std::string str;
    std::cin >> str;
    int n = str.length();
    char *s = new char[n + 2]{0};
    for (int i = 0; i < n; ++i) s[i] = str[i];
    solve(s, n);
    delete []s;
    return 0;
}
