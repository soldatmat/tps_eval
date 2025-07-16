using JSON

function prepare_input_json(sequence_id, sequence, save_path)
    data = Dict([
        "name" => sequence_id,
        "modelSeeds" => [42],
        "sequences" => [
            "protein" => Dict([
                "id" => ["A"],
                "sequence" => sequence,
            ]),
        ],
        "dialect" => "alphafold3",
        "version" => 2,
    ])

    open(save_path, "w") do f
        JSON.print(f, data, 4)
    end
end
