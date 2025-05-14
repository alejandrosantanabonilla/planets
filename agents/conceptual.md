graph TD
    subgraph User Interaction
        A[User] -- "Natural Language Query" --> B(User Query Understanding Agent / LLM Agent)
    end

    subgraph Agent Orchestration
        B -- "Identified Molecule & Property" --> C(Molecular Representation Agent)
        C -- "Molecular Representation (e.g., SMILES)" --> D(Parameter Formatting Agent)
        D -- "Formatted Parameters" --> E(GNN Execution Agent)
        E -- "Predicted Property Value" --> F(Response Generation Agent)
        F -- "Natural Language Response" --> A
    end

    subgraph Pre-trained GNN
        E -- "Molecular Representation & Property Type" --> G[Pre-trained Graph Neural Network]
        G -- "Predicted Property" --> E
    end

    subgraph External Tools & Data
        C -- "Molecule Name/Description" --> H(Molecular Databases / Cheminformatics Libraries)
        H -- "Molecular Representation" --> C
    end

    style A fill:#f9f,stroke:#333,stroke-width:2px
    style G fill:#ccf,stroke:#333,stroke-width:2px

    direction LR
