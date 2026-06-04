#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#define field 65521

uint8_t n = 11; // a, b, ~a, ~b, 1, a+b, a-b, a^b, a&b, a|b, a*b
uint8_t m = 200; // abritrary number, should be enough

#define just_a 0
#define just_b 1
#define xor 2
#define and 3
#define or 4
#define mul 5
#define add 6
#define sub 7
#define one 8
#define not_a 9
#define not_b 10

// tested, works
__attribute__((always_inline)) inline uint16_t powmod(uint64_t a, uint64_t b) {
    uint64_t acum = 1;
    while (true) {
        if (b==0) {
            return acum;
        }
        acum = (acum * a) % field; // thank you cryptohack for teaching me this
        b-=1;
    }
}

// supposed to be a helper function but makes me feel like pirate software
 __attribute__((always_inline)) inline uint16_t gf_math(uint8_t op, uint16_t a, uint16_t b) {
     uint64_t tmp = 0;
     switch (op) {
        case just_a:
            tmp = (uint64_t)a;
            break;
        case just_b:
            tmp = (uint64_t)b;
            break;
        case one:
            tmp = 1;
            break;
        case add:
            tmp = (uint64_t)a + (uint64_t)b;
            break;
        case sub:
            tmp = (uint64_t)a + field - (uint64_t)b;
            break;
        case xor:
            tmp = (uint64_t)a ^ (uint64_t)b;
            break;
        case and:
            tmp = (uint64_t)a & (uint64_t)b;
            break;
        case or:
            tmp = (uint64_t)a | (uint64_t)b;
            break;
        case mul:
            tmp = (uint64_t)a * (uint64_t)b;
            break;
        case not_a:
            tmp = (uint64_t)~a;
            break;
        case not_b:
            tmp = (uint64_t)~b;
            break;
    }
     tmp %= field;
     return (uint16_t)tmp;
}

uint16_t *generate_identity_matrix() { // wrong name?
    uint16_t* mat = malloc(n*m*sizeof(uint16_t)); 
    
    for (uint64_t a = 0; a < m; a++) {
        uint16_t tmp_a = rand();
        uint16_t tmp_b = rand();
        for (uint8_t b = 0; b < n; b++) {
            mat[a*n+b] = gf_math(b, tmp_a, tmp_b);
        }
    }

    return mat;
}

uint16_t *row_echelon_form(uint16_t *t) {
    uint16_t *ptr = (uint16_t *)malloc(n*m*sizeof(uint16_t));
    memcpy(ptr, t, n*m*sizeof(uint16_t));
    uint8_t swap_w = 0;
    uint32_t tmp_index;
    for (uint32_t a = 0; a < n; a++) {
        tmp_index = 0;
        bool found_v = false;
        for (uint32_t i = swap_w; i < m; i++) { // find non-zero coordinate in desired column
            if (ptr[n*i+a] != 0) {
                tmp_index = i;
                found_v = true;
                break;
            }
        }
        if (!found_v) {
            continue;
        }
        for (uint8_t c = 0; c < n; c++) { // swap first non-zero and the rref part
            uint16_t temp_coord = 0;
            temp_coord = ptr[swap_w*n+c];
            ptr[swap_w*n+c] = ptr[tmp_index*n+c];
            ptr[tmp_index*n+c] = temp_coord;
        }
        for (uint16_t b = (swap_w+1); b < m; b++) { // subtract from each vect
            if (ptr[b*n+a] == 0) {
                continue;
            }
            uint16_t mul_by = gf_math(mul, (ptr[b*n+a] % field), (powmod(ptr[swap_w*n+a], field-2)) % field); // modular inverse or something
            for (uint8_t c = 0; c < n; c++) {
                ptr[b*n+c] = gf_math(sub, ptr[b*n+c], gf_math(mul, mul_by, ptr[swap_w*n+c])); // subtract from full vect
            }
        }
        swap_w++; // a!=swap_w cause of 0 colums (free variables and stuff)
    }
    return ptr;
}

uint16_t *gaussian_eliminate(uint16_t *mat, uint16_t *t_v, uint16_t *pivot_rows, uint8_t pivot_tail) {
    uint16_t *s_vect = malloc(sizeof(uint16_t)*n);
    memset(s_vect, 0, sizeof(uint16_t)*n);
    for (uint8_t a = 1; a < (pivot_tail+1); a++) {
        uint16_t solve_for = 0;
        for (uint16_t i = 0; i < m; i++) {
            if (mat[i*n+pivot_rows[pivot_tail-a]] != 0) {
                solve_for = i;
            }
        }
        uint16_t mul_by = gf_math(mul, (t_v[solve_for] % field), (powmod(mat[solve_for*n+pivot_rows[pivot_tail-a]], field-2)) % field);
        
        s_vect[pivot_rows[pivot_tail-a]] = mul_by;
        for (uint16_t b = 0; b < m; b++) {
            t_v[b] = gf_math(sub, t_v[b], gf_math(mul, mat[n*b+pivot_rows[pivot_tail-a]], mul_by));
        }
    }
    return s_vect;
}

uint16_t *mat_mul(uint16_t *mat, uint16_t *vect) {
    uint16_t *r_vect = malloc(m*sizeof(uint16_t));
    uint16_t temporary = 0;
    for (uint16_t a = 0; a < m; a++) {
        temporary = 0;
        for (uint16_t b = 0; b < n; b++) {
            temporary = gf_math(add, temporary, gf_math(mul, mat[a*n+b], vect[b]));
        }
        r_vect[a] = temporary;
    }
    return r_vect;
}

uint16_t *nullspace_gen_and_sample () {
    uint16_t* matrix = generate_identity_matrix();
    uint16_t* rref_matrix = row_echelon_form(matrix);
    uint8_t cray = 0;
    uint16_t* free_cols = malloc(n*sizeof(uint16_t));
    uint8_t free_tail = 0;
    uint16_t* pivot_cols = malloc(n*sizeof(uint16_t));
    uint8_t pivot_tail = 0;

    for (uint8_t a = 0; a < n; a++) {
        for (int i = cray; i < n; i++) {
            if (rref_matrix[a*n+i] != 0) {
                cray++;
                break;
            }
            free_cols[free_tail++] = i;
            cray++;
        }
    }

    for (uint8_t a = 0; a < n; a++) {
        bool did_occur = false;
        for (int i = 0; i < free_tail; i++) {
            if (free_cols[i] == a) {
                did_occur = true;
                break;
            }
        }
        if (!did_occur) {
            pivot_cols[pivot_tail++] = a;
            did_occur = false;
        }
    }

    uint16_t *target_vect = malloc(sizeof(uint16_t)*n);

    for (uint8_t i = 0; i < n; i++) {
        for (uint8_t a = 0; a < free_tail; a++) {
            if (i == free_cols[a]) {
                target_vect[i] = rand() % field;
                break;
            } else {
                target_vect[i] = 0;
            }
        }
    }
    uint16_t *ignore_me_vect = mat_mul(rref_matrix, target_vect); // so it can be freed
    uint16_t *null_vect = gaussian_eliminate(rref_matrix, ignore_me_vect, pivot_cols, pivot_tail);
    
    for (uint8_t i = 0; i < n; i++) {
        null_vect[i] = (gf_math(sub, target_vect[i], null_vect[i]));
    }
    putchar('\n');

    free(ignore_me_vect);
    free(free_cols);
    free(pivot_cols);
    free(target_vect);
    free(matrix);
    free(rref_matrix);

    return null_vect;
}

bool null_vect_valid(uint16_t *vect) {
    for (int i = 0; i < 1000; i++) {
        uint16_t acum = 0;
        uint16_t a = rand()%field;
        uint16_t b = rand()%field;
        for (int x = 0; x < n; x++) {
            acum = gf_math(add, acum, gf_math(mul, vect[x], gf_math(x, a, b)));
        }
        if (acum != 0) {
            printf("Oops\n");
            return false;
        }
    }
    return true;
}


// try to randomly select the first term of something, then rewrite and reallocate
char *rewrite(char* input, uint64_t size) {
    try_again:
    uint64_t ptr = rand() % size;
    bool need_inc = false;
    while (input[ptr] >= '0' && input[ptr] <= '9') {
        if (ptr == 0) {
            need_inc = false;
            break;
        }
        ptr--;
        need_inc = true;
    }
    if (need_inc) {
        ptr++;
    }
    if (input[ptr] >= '0' && input[ptr] <= '9' ||
        input[ptr] >= 'a' && input[ptr] <= 'z' ||
        input[ptr] == '(') {
        uint64_t begin = ptr;
        uint64_t begin_one = begin;
        if (begin_one > 0 && input[begin_one - 1] == '*') {
            goto try_again;
        }
        int depth = 0;
        while (input[ptr] >= '0' && input[ptr] <= '9' ||
        input[ptr] >= 'a' && input[ptr] <= 'z' ||
        (input[ptr] == '(' || depth != 0)) {
            if (input[ptr] == '(') {
                depth++;
            } else if (input[ptr] == ')') {
                depth--;
            }
            if (ptr >= size) {
                goto try_again;
            }
            ptr++;
        }
        uint64_t size_var_one = ptr-begin;
        uint8_t op = input[ptr++]; // should be garaunteed to be the op
        uint8_t tmp = 0;
        switch (op) {
           case '+':
               tmp = add;
               break;
           case '-':
               tmp = sub;
               break;
           case '^':
               tmp = xor;
               break;
           case '&':
               tmp = and;
               break;
           case '|':
               tmp = or;
               break;
           default:
               goto try_again;
        }
        op = tmp;
        uint64_t begin_two = 0;
        if (input[ptr] >= '0' && input[ptr] <= '9' ||
            input[ptr] >= 'a' && input[ptr] <= 'z' ||
            input[ptr] == '(') {
            begin_two = ptr;
            int depth = 0;
            while (input[ptr] >= '0' && input[ptr] <= '9' ||
            input[ptr] >= 'a' && input[ptr] <= 'z' ||
            (input[ptr] == '(' || depth != 0)) {
                if (input[ptr] == '(') {
                    depth++;
                } else if (input[ptr] == ')') {
                    depth--;
                }
                if (ptr == size) {
                    goto try_again;
                }
                ptr++;
            }
        } else {
            goto try_again;
        }
        if (input[ptr] == '*') {
            goto try_again;
        }
        uint64_t size_var_two = ptr-(begin_two);
        uint16_t *null_vect = nullspace_gen_and_sample();
        while (!null_vect_valid(null_vect)) {
            free(null_vect);
            null_vect = nullspace_gen_and_sample();
        }



        null_vect[op] = gf_math(add, null_vect[op], 1);
        char *tmp_buff = malloc((n * (size_var_one + size_var_two + 20) + 1)*sizeof(char));
        char *var_one = strndup(input + begin_one, size_var_one); 
        char *var_two = strndup(input + begin_two, size_var_two);
        
        uint64_t tmp_buff_ptr = 0;
        
        for (uint8_t i = 0; i < n; i++) {
            if (null_vect[i] == 0) continue;
            if (tmp_buff_ptr > 0) {
                tmp_buff[tmp_buff_ptr++] = '+';
            }
            char op_chars[] = {
                [just_a] = 0, [just_b] = 0,
                [xor] = '^', [and] = '&', [or] = '|',
                [mul] = '*', [add] = '+', [sub] = '-'
            };
            if (i == just_a) {
                tmp_buff_ptr += sprintf(tmp_buff + tmp_buff_ptr,
                    "%hu*(%s)", null_vect[i], var_one);
            } else if (i == not_a) {
                tmp_buff_ptr += sprintf(tmp_buff + tmp_buff_ptr,
                    "%hu*(~%s)", null_vect[i], var_one);
            } else if (i == just_b) {
                tmp_buff_ptr += sprintf(tmp_buff + tmp_buff_ptr,
                    "%hu*(%s)", null_vect[i], var_two);
            } else if (i == not_b) {
                tmp_buff_ptr += sprintf(tmp_buff + tmp_buff_ptr,
                    "%hu*(~%s)", null_vect[i], var_two); 
            } else if (i == one) {
                tmp_buff_ptr += sprintf(tmp_buff + tmp_buff_ptr,
                    "%hu", null_vect[i]);
            } else if (op_chars[i]) {
                tmp_buff_ptr += sprintf(tmp_buff + tmp_buff_ptr,
                    "%hu*(%s%c%s)", null_vect[i], var_one, op_chars[i], var_two);
            }
        }
        tmp_buff[tmp_buff_ptr] = '\0';
        
        uint64_t orig_len = strlen(input);
        uint64_t seg_start = begin_one;           
        uint64_t seg_end = ptr;                 
        uint64_t seg_len = seg_end - seg_start;
        
        char *new_input = malloc(orig_len - seg_len + tmp_buff_ptr + 3); 
        uint64_t w = 0;
        memcpy(new_input + w, input, seg_start);
        w += seg_start;
        new_input[w++] = '(';
        memcpy(new_input + w, tmp_buff, tmp_buff_ptr);
        w += tmp_buff_ptr;
        new_input[w++] = ')';
        memcpy(new_input + w, input + seg_end, orig_len - seg_end + 1);
        
        free(var_one);
        free(var_two);
        free(tmp_buff);
        free(null_vect);
        free(input);
        
        return new_input;
    } else {
        goto try_again;
    }
}

    

void main () {
    char *input = 0;
    size_t cap = 0;
    printf("Equation: ");
    uint64_t size = getline(&input, &cap, stdin);
    printf("Effort (int): ");
    uint32_t effort = 0;
    scanf("%u", &effort);
    printf("Seed (zero for random): ");
    uint32_t seed = 0;
    scanf("%u", &seed);
    if (seed == 0) {
        seed = rand();
        printf("Seed: %u\n", seed);
    }
    srand(seed);
    while (effort--) {
        input = rewrite(input, strlen(input));
    }
    printf("%s", input);

}

