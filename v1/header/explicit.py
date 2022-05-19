
floats = ["float", "double", "long double"]
debug = ["true", "false"]
par = ["qsl::seq", "qsl::omp", "qsl::opt"]
sims = ["qsl::basic", "qsl::resize", "qsl::number"]

for s in sims:
    for f in floats:
        for d in debug:
            for p in par:
                print(f"template class {s}<{f}, {d}, {p}>")


