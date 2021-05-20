for file in dags/*/*.dag
do
    echo "$(basename ${file}):"
    ./c3po "${file}"
    echo "---"
done
