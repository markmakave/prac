#include <stdio.h>
#include <stdlib.h>
#include <papi.h>

struct item
{
    double val;
    struct item *next;
};

struct item *
gen_list(int n)
{
    struct item *res = NULL;
    int i;
    for (i = 0; i < n; ++i) {
        struct item *new_el;

        new_el = calloc(1, sizeof(*new_el));
        if (!new_el) {
            fputs("memory allocation error\n", stderr);
            return NULL;
        }
        new_el->next = res;
        new_el->val = (double) rand();
        res = new_el;
    }
    return res;
}

void
list_free(struct item *l)
{
    struct item *p = l;
    while(p) {
        struct item *tmp = p;
        p = p->next;
        free(tmp);
    }
}

double
list_sum(struct item *l)
{
    struct item *p = l;
    double sum = 0.0;

    while(p) {
        sum += p->val;
        p = p->next;
    }
    return sum;
}


void
gen_array(double *a, int size)
{
    int i;
    for (i = 0; i < size; ++i) {
        a[i] = (double) rand();
    }
}

double
array_sum(double *a, int size)
{
    double res = 0.0;
    int i;
    for (i = 0; i < size; ++i) {
        res += a[i];
    }
    return res;
}

void
code1(int n)
{
    
    double *a = NULL, sum;

    a = calloc(n, sizeof(a[0]));
    if (!a) {
        fputs("failed to allocate array a\n", stderr);
        return;
    }

    srand(323);

    gen_array(a, n);
    sum = array_sum(a, n);

    printf("code1: sum = %f\n", sum);

    free(a);
}

void
code2(int n)
{
    struct item *l = NULL;
    double sum;
    
    srand(323);

    l = gen_list(n);

    if (!l) {
        fputs("failed to generate list\n", stderr);
        return;
    }

    sum = list_sum(l);

    printf("code2: sum = %f\n", sum);

    list_free(l);
}


int
main(int argc, char **argv)
{
    int retval;
    int event_set = PAPI_NULL;
    long long cm1, cm2;
    
    int n = 1024 * 1024;

    retval = PAPI_library_init(PAPI_VER_CURRENT);
    if (retval != PAPI_VER_CURRENT) {
        fputs("failed to initialize PAPI\n", stderr);
        return 1;
    }

    if (PAPI_create_eventset(&event_set) != PAPI_OK) {
        fputs("unable to create event set\n", stderr);
        return 2;
    }

    if (PAPI_add_event(event_set, PAPI_L1_TCM) != PAPI_OK) {
        fputs("can't add event\n", stderr);
        return 3;
    }

    if (PAPI_start(event_set) != PAPI_OK) {
        fputs("failed to start counters\n", stderr);
        return 4;
    }

    code1(n);

    if (PAPI_stop(event_set, &cm1) != PAPI_OK) {
        fputs("error in stop counters 1\n", stderr);
        return 5;
    }

    if (PAPI_start(event_set) != PAPI_OK) {
        fputs("failed to start counters\n", stderr);
        return 6;
    }

    code2(n);

    if (PAPI_stop(event_set, &cm2) != PAPI_OK) {
        fputs("error in stop counters 1\n", stderr);
        return 7;
    }

    printf("cache misses for code 1 = %lld\n"
           "cache misses for code 2 = %lld\n",
           cm1, cm2);
    
    return 0;
}
