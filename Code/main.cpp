#include "stdafx.h"

int main(int argc, char * argv[]) {
    const TArgument Arg(argc, argv);
    string adjpath, corepath, pagepath;
    if (Arg._graphname == "dblp") {
        adjpath = "spreprocess/com-dblp-adj.file";
        corepath = "spreprocess/com-dblp-core.file";
        pagepath = "spreprocess/com-dblp-page.file";
    }
    else if (Arg._graphname == "dblp_case") {
        adjpath = "spreprocess/dblp-case-adj.file";
        corepath = "spreprocess/dblp-case-core.file";
        pagepath = "spreprocess/dblp-case-page.file";
    }
    else if (Arg._graphname == "domain_case") {
        adjpath = "spreprocess/domain-case-adj.file";
        corepath = "spreprocess/domain-case-core.file";
        pagepath = "spreprocess/domain-case-page.file";
    }
    else if (Arg._graphname == "domain_pub") {
        adjpath = "spreprocess/domain-pub-adj.file";
        corepath = "spreprocess/domain-pub-core.file";
        pagepath = "spreprocess/domain-pub-page.file";
    }
    else if (Arg._graphname == "domain_label") {
        adjpath = "spreprocess/domain-label-adj.file";
        corepath = "spreprocess/domain-label-core.file";
        pagepath = "spreprocess/domain-label-page.file";
    }
    else if (Arg._graphname == "topic_case") {
        adjpath = "spreprocess/topic-case-adj.file";
        corepath = "spreprocess/topic-case-core.file";
        pagepath = "spreprocess/topic-case-page.file";
    }
    else if (Arg._graphname == "topic_pub") {
        adjpath = "spreprocess/topic-pub-adj.file";
        corepath = "spreprocess/topic-pub-core.file";
        pagepath = "spreprocess/topic-pub-page.file";
    }
    else if (Arg._graphname == "topic_label") {
        adjpath = "spreprocess/topic-label-adj.file";
        corepath = "spreprocess/topic-label-core.file";
        pagepath = "spreprocess/topic-label-page.file";
    }
    else if (Arg._graphname == "orkut") {
        adjpath = "spreprocess/com-orkut-adj.file";
        corepath = "spreprocess/com-orkut-core.file";
        pagepath = "spreprocess/com-orkut-page.file";
    }
    else if (Arg._graphname == "email") {
        adjpath = "spreprocess/email-adj.file";
        corepath = "spreprocess/email-core.file";
        pagepath = "spreprocess/email-page.file";
    }
    else if (Arg._graphname == "youtube") {
        adjpath = "spreprocess/com-youtube-adj.file";
        corepath = "spreprocess/com-youtube-core.file";
        pagepath = "spreprocess/com-youtube-page.file";
    }
    else if (Arg._graphname == "lj") {
        adjpath = "spreprocess/com-lj-adj.file";
        corepath = "spreprocess/com-lj-core.file";
        pagepath = "spreprocess/com-lj-page.file";
    }
    else if (Arg._graphname == "fs") {
        adjpath = "spreprocess/com-friendster-adj.file";
        corepath = "spreprocess/com-friendster-core.file";
        pagepath = "spreprocess/com-friendster-page.file";
    }
//    string path = "D:\\Git\\Weight_Community\\preprocess\\com-dblp.txt";

    if (Arg._mode == "top") {
        TopRComm * app = new TopRComm();
        if (Arg._func == "max") {
            if (Arg._algName == "naive") {
                app->naive_max_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "improve") {
                app->cons_max_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._p);
            }
        }
        else if (Arg._func == "min") {
            if (Arg._algName == "naive") {
                app->naive_min_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "improve") {
                app->cons_min_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._p);
            }
            else if (Arg._algName == "non") {
                app->non_min_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._p);
            }
        }
        else if (Arg._func == "free_sum") {
            if (Arg._algName == "naive") {
                app->naive_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "improve") {
                app->improved_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "approx") {
                app->approx_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._eps);
            }
        }
        else if (Arg._func == "free_avg") {
            if (Arg._algName == "random") {
                app->improved_random_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "improve") {
                app->improved_greedy_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "climb") {
                app->improved_climb_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "double") {
                app->improved_double_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
            else if (Arg._algName == "cut") {
                app->cut_double_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r);
            }
        }
        else if (Arg._func == "sum") {
            if (Arg._algName == "naive") {
                app->naive_cons_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r);
            }
            else if (Arg._algName == "improve") {
                app->improved_cons_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r);
            }
            else if (Arg._algName == "approx") {
                app->approx_cons_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._eps);
            }
        }
        else if (Arg._func == "non_sum") {
            if (Arg._algName == "improve") {
                app->non_cons_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
        }
        else if (Arg._func == "size_sum") {
            if (Arg._algName == "climb") {
                app->climb_size_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._s);
            }
            else if (Arg._algName == "degree") {
                app->degree_size_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._s);
            }
            else if (Arg._algName == "casual") {
                app->casual_size_sum_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._s);
            }
        }
        else if (Arg._func == "avg") {
            if (Arg._algName == "random") {
                app->improved_cons_random_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
            else if (Arg._algName == "improve") {
                app->improved_cons_greedy_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
            else if (Arg._algName == "climb") {
                app->improved_cons_climb_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
            else if (Arg._algName == "double") {
                app->improved_cons_double_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
            else if (Arg._algName == "cut") {
                app->cut_cons_double_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
            else if (Arg._algName == "exact") {
                app->improved_cons_exact_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath,
                        Arg._k, Arg._r, Arg._p);
            }
        }
        else if (Arg._func == "size_avg") {
            if (Arg._algName == "climb") {
                app->climb_size_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._s);
            }
            else if (Arg._algName == "degree") {
                app->degree_size_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._s);
            }
            else if (Arg._algName == "casual") {
                app->casual_size_avg_global_topr(Arg._graphname, adjpath, corepath, pagepath, Arg._k, Arg._r, Arg._s);
            }
        }
    }

    std::cout << "Hello3, World!" << std::endl;
    return 0;
}
