class WorkflowChecker:
    def __init__(self):
        self.steps = {
            'fetch_term_names': False,
            'fetch_term_products': False,
            'create_products': False,
            'fetch_uniprotid': False,
            'prune_products': False,
            'fetch_uniprot_infos': False,
            'score_products': False,
            'fetch_mRNA': False,
            'predict_miRNA': False,
            'score_miRNA': False,
        }
    
    def mark_step_done(self, step_name):
        if step_name in self.steps:
            self.steps[step_name] = True
        else:
            raise ValueError(f"Invalid step name: {step_name}")
    
    def is_step_done(self, step_name):
        if step_name in self.steps:
            return self.steps[step_name]
        else:
            raise ValueError(f"Invalid step name: {step_name}")