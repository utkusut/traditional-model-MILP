#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

int main(int argc, char** argv) {
    IloEnv env;
    try {
        // Define the model
        IloModel model(env);

        // Parameters
        const int numBlocks = 28;
        const int numPeriods = 5;
        // Block values
        IloNumArray blockValues(env, 28, -100, -100, 200, -100, -100, -100, -100, -100, 1000, 400, 600, 700, 300, -100, -100, -100, 700, 800, -100, -100, -100, -100, -100, -100, 1000, -100, -100, -100);
        IloNumArray weights(env, numBlocks, 1000, 1000, 2000, 2000, 1000, 1000, 1000, 1000, 2500, 2500, 2500, 2500, 2500, 1000, 1000, 1000, 3000, 3000, 3000, 1000, 1000, 1000, 1000, 1000, 3500, 1000, 1000, 1000); // Fill in the block weights
        IloNumArray miningCapacity(env, numPeriods, 8000, 8000, 8000, 8000, 8000);
        IloNumArray processingCapacity(env, numPeriods, 8000, 8000, 8000, 8000, 8000);

        // Decision variables
        IloArray<IloBoolVarArray> X(env, numBlocks);
        for (int i = 0; i < numBlocks; i++) {
            X[i] = IloBoolVarArray(env, numPeriods);
        }

        IloExpr obj(env);
        for (int b = 0; b < numBlocks; b++) {
            for (int t = 0; t < numPeriods; t++) {
                obj += (blockValues[b] / pow(1 + discountRate, t + 1)) * X[b][t];
             }
}
model.add(IloMaximize(env, obj));

        // Constraints

       // Reserve constraint
    for (int b = 0; b < numBlocks; b++) {
        IloExpr reserveExpr(env);
        for (int t = 0; t < numPeriods; t++) {
            reserveExpr += X[b][t];
            }
                    model.add(reserveExpr <= 1);
            }


         // Mining capacity constraint
        for (int t = 0; t < numPeriods; t++) {
            IloExpr miningExpr(env);
            for (int b = 0; b < numBlocks; b++) {
                miningExpr += weights[b] * X[b][t];
            }
            model.add(miningExpr <= miningCapacity[t]);
        }

        // Processing capacity constraint
    for (int t = 0; t < numPeriods; t++) {
     IloExpr processingExpr(env);
     for (int b = 0; b < numBlocks; b++) {
           if (blockValues[b] > 0) {
              processingExpr += weights[b] * X[b][t];
          }
        }
    model.add(processingExpr <= processingCapacity[t]);
        }


        // Precedence constraints
        IloArray<IloNumArray> precedence(env, 28);

            precedence[7] = IloNumArray(env, 2, 0, 1);     // Block 8 has predecessors 1 and 2
            precedence[8] = IloNumArray(env, 3, 0, 1, 2);  // Block 9 has predecessors 1, 2, and 3
            precedence[9] = IloNumArray(env, 3, 1, 2, 3);  // Block 10 has predecessors 2, 3, and 4
            precedence[10] = IloNumArray(env, 3, 2, 3, 4); // Block 11 has predecessors 3, 4, and 5
            precedence[11] = IloNumArray(env, 3, 3, 4, 5); // Block 12 has predecessors 4, 5, and 6
            precedence[12] = IloNumArray(env, 3, 4, 5, 6); // Block 13 has predecessors 5, 6, and 7
            precedence[13] = IloNumArray(env, 2, 5, 6);    // Block 14 has predecessors 6 and 7
            precedence[14] = IloNumArray(env, 2, 7, 8);    // Block 15 has predecessors 8 and 9
            precedence[15] = IloNumArray(env, 3, 7, 8, 9); // Block 16 has predecessors 8, 9, and 10
            precedence[16] = IloNumArray(env, 3, 8, 9, 10);// Block 17 has predecessors 9, 10, and 11
            precedence[17] = IloNumArray(env, 3, 9, 10, 11);// Block 18 has predecessors 10, 11, and 12
            precedence[18] = IloNumArray(env, 3, 10, 11, 12);// Block 19 has predecessors 11, 12, and 13
            precedence[19] = IloNumArray(env, 3, 11, 12, 13);// Block 20 has predecessors 12, 13, and 14
            precedence[20] = IloNumArray(env, 2, 12, 13);  // Block 21 has predecessors 13 and 14
            precedence[21] = IloNumArray(env, 2, 14, 15);  // Block 22 has predecessors 15 and 16
            precedence[22] = IloNumArray(env, 3, 14, 15, 16); // Block 23 has predecessors 15, 16, and 17
            precedence[23] = IloNumArray(env, 3, 15, 16, 17); // Block 24 has predecessors 16, 17, and 18
            precedence[24] = IloNumArray(env, 3, 16, 17, 18); // Block 25 has predecessors 17, 18, and 19
            precedence[25] = IloNumArray(env, 3, 17, 18, 19); // Block 26 has predecessors 18, 19, and 20
            precedence[26] = IloNumArray(env, 3, 18, 19, 20); // Block 27 has predecessors 19, 20, and 21
            precedence[27] = IloNumArray(env, 2, 19, 20);  // Block 28 has predecessors 20 and 21

        for (int b = 0; b < 28; b++) {
    int n = precedence[b].getSize(); // Number of predecessor blocks for block b
    for (int t = 0; t < numPeriods; t++) {
        IloExpr expr(env);
        for (int i = 0; i < n; i++) {
            int b_prev = precedence[b][i]; // Predecessor block
            for (int t_prev = 0; t_prev <= t; t_prev++) {
                expr += X[b_prev][t_prev];
            }
        }
        model.add(n * X[b][t] - expr <= 0);
    }
}


// Solve the model
IloCplex cplex(model);
if (cplex.solve()) {
    env.out() << "Solution status: " << cplex.getStatus() << endl;
    env.out() << "Objective value: " << cplex.getObjValue() << endl;

    // Print the values of X to see which blocks are mined in which years
    for (int b = 0; b < 28; b++) {
        for (int t = 0; t < numPeriods; t++) {
            if (cplex.getValue(X[b][t]) > 0.5) { // Assuming binary decision variables
                env.out() << "Block " << (b + 1) << " should be mined in year " << (t + 1) << endl;
            }
        }
    }

    // If you have capacity constraints and want to see the used capacities, you can similarly retrieve the values of the associated expressions.

} else {
    env.out() << "Failed to optimize." << endl;
}

} catch (IloException& e) {
    cerr << "Exception caught: " << e << endl;
} catch (...) {
    cerr << "Unknown exception caught" << endl;
}

env.end();
return 0;
