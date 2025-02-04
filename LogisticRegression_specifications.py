from q2_feature_classifier.classifier import spec_from_pipeline
from sklearn.feature_extraction.text import HashingVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from json import dumps, dump

steps = [
    ('feat_ext',
     HashingVectorizer(
         analyzer='char_wb',
         # increased from the values used in the forum, because apparently
         # a higher number is better for DNA/RNA sequences
         n_features=65536,
         ngram_range=[8,8],
         alternate_sign=False)),
    ('classify',
     LogisticRegression(max_iter=1000))
]

pipeline = Pipeline(steps=steps)
spec = spec_from_pipeline(pipeline)
classifier_specification = dumps(spec, indent=2)
print(classifier_specification)

with open('LogReg_classifier_specification.json', 'w') as file:
    dump(spec, file, indent=2)
