#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <minhash.h>
#include <configuration.h>

struct minhash_configuration conf = {
    .sketch_size = 128,          /// Number of hash functions / sketch size
    .prime_modulus = (1ULL << 31) - 1,       /// Large prime for hashing (M)
    .hash_type = 0,        /// ID for hash function pointer
    .init_size = 0,                 /// Initial elements to insert (optional)
    .k = 5,
    .N = 0,
    .b = 0,
};

int main(int argc, const char* argv[]) {

    if (argc < 6) {
        fprintf(stderr, 
            "Parameter error! Make sure to pass <N_total: (int) total unique element budget, sketch_size: (int) size of the sketch, init_size: (int) starting size of set, l_percent: (double) percent to omit from B, true_jaccard: (double) targeted true jaccard>\n");
        exit(1);
    }

    char *endptr;
    long n_total = strtol(argv[1], &endptr, 10);
    if (*endptr != '\0' || n_total <= 0){
        fprintf(stderr, "n_total must be greater than zero!\n");
        exit(1);
    }

    long ssize = strtol(argv[2], &endptr, 10);
    if (*endptr != '\0' || ssize <= 0) {
        fprintf(stderr, "Size of sketch must be greater than zero!\n");
        exit(1);
    }

    long startsize = strtol(argv[3], &endptr, 10);
    if (*endptr != '\0' || startsize <= 0) {
        fprintf(stderr, "Starting size of set must be greater than zero!\n");
        exit(1);
    }
    
    double l_percent = strtod(argv[4], &endptr);
    if (*endptr != '\0' || l_percent < 0.0 || l_percent > 1.0) {
        fprintf(stderr, "l_percent must be a valid float between 0.0 and 1.0!\n");
        exit(1);
    }

    double true_jaccard = strtod(argv[5], &endptr);
    if (*endptr != '\0' || true_jaccard < 0.0 || true_jaccard > 1.0) {
        fprintf(stderr, "true_jaccard must be a valid float between 0.0 and 1.0!\n");
        exit(1);
    }
    
    // --- Jaccard Math Set Boundaries ---
    // J = I / (I + 2U) -> Solving for I and U given n_total = I + 2U
    long intersection_size = (long)round(n_total * (true_jaccard / (1.0 + true_jaccard)));
    long unique_per_set = (n_total - intersection_size) / 2;
    
    // Recalculate actual true Jaccard due to integer rounding constraints
    double actual_target_jaccard = (double)intersection_size / (intersection_size + 2 * unique_per_set);

    printf("Targeted Jaccard: %f | Actual Set Jaccard (due to rounding): %f\n", true_jaccard, actual_target_jaccard);
    printf("Elements Allocation -> Intersection: %ld, Unique Target: %ld, Unique A: %ld\n", 
           intersection_size, unique_per_set, unique_per_set);
           
    conf.sketch_size = (uint64_t) ssize;
    if (startsize > 0) conf.init_size = (uint64_t) startsize;

    read_configuration(conf);

    minhash_sketch *target_sketch;
    minhash_sketch *sketch_A;
    minhash_sketch *sketch_B;

    void *hash_functions = hash_functions_init(conf.hash_type, conf.sketch_size, conf.prime_modulus, conf.k);

    minhash_init(&target_sketch, hash_functions, conf.sketch_size, conf.init_size, conf.hash_type);
    minhash_init(&sketch_A, hash_functions, conf.sketch_size, conf.init_size, conf.hash_type);
    minhash_init(&sketch_B, hash_functions, conf.sketch_size, conf.init_size, conf.hash_type);

    long element_pool = startsize;

    // 1. Insert Shared Intersection Elements into all sketches
    // For Sketch B, we track its total planned inserts so we can drop the last l_percent
    long total_planned_A_inserts = intersection_size + unique_per_set;
    long b_cutoff = (long)(total_planned_A_inserts * (1.0 - l_percent));
    long b_inserted_count = 0;

    long i;
    for (i = 0; i < intersection_size; i++) {
        insert(target_sketch, element_pool);
        insert(sketch_A, element_pool);
        
        if (b_inserted_count < b_cutoff) {
            insert(sketch_B, element_pool);
            b_inserted_count++;
        }
        element_pool++;
    }

    // 2. Insert Unique elements into Target Sketch only
    for (i = 0; i < unique_per_set; i++) {
        insert(target_sketch, element_pool);
        element_pool++;
    }

    // 3. Insert Unique elements into Sketch A (and Sketch B, subject to cutoff)
    for (i = 0; i < unique_per_set; i++) {
        insert(sketch_A, element_pool);
        
        if (b_inserted_count < b_cutoff) {
            insert(sketch_B, element_pool);
            b_inserted_count++;
        }
        element_pool++;
    }

    // Query similarities against Target sketch
    float est_A = query(target_sketch, sketch_A);
    float est_B = query(target_sketch, sketch_B);

    // Compute absolute errors relative to the actual generated set Jaccard
    float error_A = fabs(est_A - actual_target_jaccard);
    float error_B = fabs(est_B - actual_target_jaccard);

    // Output results
    printf("\n--- Similarity Evaluation Results ---\n");
    printf("Actual Set Jaccard: %f\n", actual_target_jaccard);
    printf("Estimated Jaccard A: %f | Absolute Error A: %f\n", est_A, error_A);
    printf("Estimated Jaccard B: %f | Absolute Error B: %f\n", est_B, error_B);

    // Securely deallocate memory
    minhash_free(target_sketch);
    minhash_free(sketch_A);
    minhash_free(sketch_B);

    return 0;
}
