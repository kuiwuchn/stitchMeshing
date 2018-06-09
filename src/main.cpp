#include "common.h"
#include "viewer.h"
#include "batch.h"
/* Force usage of discrete GPU on laptops */
NANOGUI_FORCE_DISCRETE_GPU();

int main(int argc, char **argv) {
    std::vector<std::string> args;
    bool fullscreen = false, help = false;
    int nprocs = -1;
    #if defined(__APPLE__)
        bool launched_from_finder = false;
    #endif

		bool Batch_Process = false;
		uint32_t  dim = 2;
		char batchIutput[300], batchOutput[300];
		Float scale = 3, tlen=1.0; uint32_t smooth_iter = 10;
    try {
        for (int i=1; i<argc; ++i) {
            if (strcmp("--fullscreen", argv[i]) == 0 || strcmp("-F", argv[i]) == 0) {
                fullscreen = true;
            } else if (strcmp("--help", argv[i]) == 0 || strcmp("-h", argv[i]) == 0) {
                help = true;
            } else if (strcmp("--threads", argv[i]) == 0 || strcmp("-t", argv[i]) == 0) {
                if (++i >= argc) {
                    cerr << "Missing thread count!" << endl;
                    return -1;
                }
                nprocs = str_to_uint32_t(argv[i]);
#if defined(__APPLE__)
            } else if (strncmp("-psn", argv[i], 4) == 0) {
                launched_from_finder = true;
#endif
            }
			else if (strcmp("-batch", argv[i]) == 0 || strcmp("-b", argv[i]) == 0) {
				Batch_Process = true;
			}
			else if (strcmp("--dimension", argv[i]) == 0 || strcmp("-d", argv[i]) == 0) {
				if (++i >= argc) {
					cerr << "Missing tet-gen density argument!" << endl;
					return -1;
				}
				dim = str_to_uint32_t(argv[i]);
			}
			else if (strcmp("--tlenratio", argv[i]) == 0 || strcmp("-tl", argv[i]) == 0) {
				if (++i >= argc) {
					cerr << "Missing tet-gen density argument!" << endl;
					return -1;
				}
				tlen = str_to_float(argv[i]);
			}
			else if (strcmp("--smooth", argv[i]) == 0 || strcmp("-S", argv[i]) == 0) {
				if (++i >= argc) {
					cerr << "Missing smoothing iteration count argument!" << endl;
					return -1;
				}
				smooth_iter = str_to_uint32_t(argv[i]);
			}
			else if (strcmp("--scale", argv[i]) == 0 || strcmp("-s", argv[i]) == 0) {
				if (++i >= argc) {
					cerr << "Missing scale argument!" << endl;
					return -1;
				}
				scale = str_to_float(argv[i]);
			}
			else if (strcmp("--iutput", argv[i]) == 0 || strcmp("-i", argv[i]) == 0) {
				if (++i >= argc) {
					cerr << "Missing batch mode output file argument!" << endl;
					return -1;
				}
				sprintf(batchIutput, "%s", argv[i]);
			}
			else if (strcmp("--output", argv[i]) == 0 || strcmp("-o", argv[i]) == 0) {
				if (++i >= argc) {
					cerr << "Missing batch mode output file argument!" << endl;
					return -1;
				}
				sprintf(batchOutput, "%s", argv[i]);
			}
			else {
                if (strncmp(argv[i], "-", 1) == 0) {
                    cerr << "Invalid argument: \"" << argv[i] << "\"!" << endl;
                    help = true;
                }
                args.push_back(argv[i]);
            }
        }
    } catch (const std::exception &e) {
        cout << "Error: " << e.what() << endl;
        help = true;
    }
	
	nprocs = -1;
    tbb::task_scheduler_init init(nprocs == -1 ? tbb::task_scheduler_init::automatic : nprocs);

	if (Batch_Process) {
		batch_process(batchIutput, batchOutput, dim, tlen, scale, smooth_iter);
		return EXIT_SUCCESS;
	}
	
    try {
        nanogui::init();

        #if defined(__APPLE__)
            if (launched_from_finder)
                nanogui::chdir_to_bundle_parent();
        #endif

		std::string filename;
        if (args.size() > 0)
            filename = args[0];
        nanogui::ref<Viewer> viewer = new Viewer(filename, fullscreen);
        viewer->setVisible(true);

        nanogui::mainloop();

        nanogui::shutdown();
    } catch (const std::runtime_error &e) {
        std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
        #if defined(_WIN32)
            MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
        #else
            std::cerr << error_msg << endl;
        #endif
        return -1;
    }

    return EXIT_SUCCESS;
}
