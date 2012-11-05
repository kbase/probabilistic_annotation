if [ `command -v dos2unix 2>/dev/null` ]; then
    dos2unix ProbabilisticAnnotation.spec
fi

base=probabilistic_annotation
compile_typespec                \
        -impl ${base}Impl           \
        -service ${base}Server      \
        -psgi ${base}.psgi          \
        -client ${base}Client       \
        -js ${base}Client           \
        -py ${base}Client           \
        ProbabilisticAnnotation.spec service