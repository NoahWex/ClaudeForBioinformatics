# Exercise 1: Project Scaffolding

**Objective**: Use Claude to create a well-organized project structure for your analysis.

**Time**: ~35 minutes

---

## Background

A consistent project structure makes analyses reproducible and shareable. In this exercise, you'll use Claude to create a standard bioinformatics project layout.

## Your Task

### Part 1: Create the directory structure (15 min)

Start a conversation with Claude and ask it to help you create a project structure. Here's an example prompt to get you started:

```
I'm starting a new scRNA-seq integration project. Please help me create a project
directory structure in ./working/ with folders for:
- scripts (R analysis scripts)
- data (input data)
- outputs (processed data, organized by date)
- reports (Quarto/Rmd reports)
- figures (publication-quality figures)

Also create a basic README.md explaining the project structure.
```

**Watch for**: Claude using the `Bash` tool to create directories and `Write` tool to create the README.

### Part 2: Create the manifest file (15 min)

Now ask Claude to help you create a manifest configuration file:

```
Please create a manifest.yaml file in ./working/ for my project. I'll be analyzing
a breast tissue scRNA-seq dataset:
- gray.rds (Gray et al., ~53K cells, 16 donors)

The data is located at:
/share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/

Include sections for:
- Project metadata (name, author, date)
- Data source with file path
- QC parameters (placeholder values for now)
- Output directories
```

**Watch for**: Claude using `Write` to create the YAML file with proper structure.

### Part 3: Verify your structure (5 min)

Ask Claude to show you what was created:

```
Please show me the current directory structure of ./working/ and display the
contents of the manifest.yaml file.
```

**Watch for**: Claude using `Bash` (ls/tree) and `Read` tools.

---

## Expected Outcome

Your `working/` directory should look like:

```
working/
├── README.md
├── manifest.yaml
├── scripts/
├── data/
├── outputs/
├── reports/
└── figures/
```

## Discussion Questions

1. What tools did Claude use to accomplish each task?
2. How would you modify the prompt if the structure wasn't quite what you wanted?
3. What's the advantage of using a manifest file over hardcoded paths?

## Hints

- Be specific about what you want, but don't over-specify
- If Claude creates something you don't like, tell it what to change
- You can reference the template in `../participant_workspace/` if needed

## Next Steps

Once your project structure is ready, proceed to [Exercise 2: Data Loading](02_data_loading.md).
