include("prepare_input.jl")

using Pkg

cd(@__DIR__)
# TODO create julia project in tps_eval
# Pkg.activate(".") # TODO uncomment
Pkg.activate("../../../../terpene_generation/src") # TODO delete
#Pkg.instantiate()

# --- Script ---------------------------------------------------------
using ArgParse

function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table s begin
        "sequence_id"
            help = "ID of the sequence"
        "sequence"
            help = "Amino acid sequence"
        "save_path"
            help = "Path to save the output JSON"
    end
    return parse_args(s)
end

args = parse_arguments()
sequence_id = args["sequence_id"]
sequence = args["sequence"]
save_path = args["save_path"]

prepare_input_json(
    sequence_id,
    sequence,
    save_path,
)
