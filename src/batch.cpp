/*
    batch.cpp -- command line interface to Instant Meshes

    This file is part of the implementation of

        Instant Field-Aligned Meshes
        Wenzel Jakob, Daniele Panozzo, Marco Tarini, and Olga Sorkine-Hornung
        In ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2015)

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE.txt file.
*/

#include "batch.h"
#include "meshio.h"
#include "hierarchy.h"
#include "optimizer.h"
#include "dedge.h"
#include "subdivide.h"
#include "timer.h"

void batch_process(char *input, char *output,
	uint32_t dimension, Float tlen, Float scale, int smooth_iter) {
	//cout << endl;
	//cout << "Running in batch mode:" << endl;
	//cout << "   Input file             = " << input << endl;
	//cout << "   Output file            = " << output << endl;
	char patho[1024];

	MultiResolutionHierarchy mRes;
	Optimizer *mOptimizer;

	Timer<> timer;
	timer.beginStage("data pre-processing");
	mRes.load(input);

	MeshStats stats = mRes.compute_mesh_stats(mRes.F(), mRes.V(0));
	mRes.diagonalLen *= 0.0025;
	if(dimension == 2 && (stats.mMaximumEdgeLength > mRes.diagonalLen)) {
	//if (dimension ==2 && ((stats.mMaximumEdgeLength * 4 / stats.mAverageEdgeLength > scale || stats.mMaximumEdgeLength > stats.mAverageEdgeLength * 1.5))) {
		cout << "Input mesh is too coarse for the desired output edge length "
			"(max input mesh edge length=" << stats.mMaximumEdgeLength
			<< "), subdividing .." << endl;
		VectorXu V2E, E2E;
		VectorXb boundary, nonManifold;
		build_dedge(mRes.F(), mRes.V(0), V2E, E2E, boundary, nonManifold);
		subdivide(mRes.F(), mRes.V(0), V2E, E2E, boundary, nonManifold, mRes.diagonalLen);

		//sprintf(patho, "%s%s", output, "_refined.obj");
		//write_surface_mesh_OBJ(mRes.V(0), mRes.F(), patho);
	}
	else return;

	//mRes.laplacian_smoothing(mRes.V(0), 3);

	//sprintf(patho, "%s%s", output, "_smoothed.obj");
	//write_surface_mesh_OBJ(mRes.V(0), mRes.F(), patho);

	mRes.build();
	
	mRes.setScale(scale);
	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	timer.beginStage("rosy optimization");

	mOptimizer = new Optimizer(mRes);
	mOptimizer->setMaxIterations(smooth_iter);
	mOptimizer->setOptimizeOrientations(true);
	mOptimizer->notify();
	mOptimizer->wait();

	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	timer.beginStage("posy optimization");

	mOptimizer->setOptimizePositions(true);
	mOptimizer->notify();
	mOptimizer->wait();

	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	mOptimizer->shutdown();

	timer.beginStage("mesh extraction");

	mRes.re_color = true;
	mRes.splitting = false;
	mRes.decomposes = true;
	mRes.doublets = true;
	mRes.triangles = true;

	if (!mRes.tetMesh()) {
		mRes.re_color = true;
		mRes.splitting = true;
		mRes.decomposes = true;
		mRes.doublets = true;
		mRes.triangles = true;

		mRes.meshExtraction2D();
	}
	else {
		mRes.re_color = true;
		mRes.splitting = true;
		mRes.meshExtraction3D();//mRes.extractTet();
	}

	timer.endStage();
	mRes.sta.timings.push_back(timer.value());

	sprintf(patho, "%s%s", output, "_surout.obj");
	write_surface_mesh_OBJ(mRes.mV_tag, mRes.F_tag, patho);
}
